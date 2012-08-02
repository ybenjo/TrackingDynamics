#include "func.hpp"

typedef std::pair<int, int> key;

class MM{
public:
  MM(){
  }

  MM(double alpha, double lambda, double k_max){
    this->alpha = alpha;
    this->lambda = lambda;
    this->k_max = k_max;
  }

  int set_word(std::string w){
    for(int i = 0; i < words.size(); ++i){
      if(words.at(i) == w){
	return i;
      }
    }
    words.push_back(w);
    return words.size();
  }

  // format
  // timestamp \t word \t word \t word \t word ...
  void set_document(const std::vector<std::string>& doc){
    std::unordered_map<int, double> counts;
    int doc_id = documents.size();
    timestamps[doc_id] = atof(doc.at(0).c_str());

    for(std::vector<std::string>::const_iterator w = doc.begin() + 1; w != doc.end(); ++w){
      int word_id = set_word(*w);
      counts[word_id]++;
    }
    documents.push_back(counts);
  }

  // update prob_doc
  double calc_prob_doc(int doc_id,
		       const std::unordered_map<int, double>& mu,
		       const std::unordered_map<int, double>& sigma
		       ){
    double sum = 0.0;

    std::unordered_map<int, double> x_minus_mu;
    std::unordered_map<int, double>::const_iterator i;
    std::unordered_map<int, double>::iterator j;

    // calc ( (x_1 - mu_1)^2, (x_2 - mu_2)^2, ... )
    for(j = documents.at(doc_id).begin(); j !=documents.at(doc_id).end(); ++j){
      x_minus_mu[(*j).first] += (*j).second;
    }
    
    for(i = mu.begin(); i != mu.end(); ++i){
      x_minus_mu[(*i).first] -= (*i).second;
    }
    
    for(j = x_minus_mu.begin(); j != x_minus_mu.end(); ++j){
      x_minus_mu[(*j).first] = pow((*j).second, 2);
    }

    for(i = sigma.begin(); i != sigma.end(); ++i){
      int key = (*i).first;
      double val = (*i).second;
      sum += val * x_minus_mu[key];
    }

    // Ignore denominator in normal distribution((2pi)^(d/2) * (\Sigma)^(1/2))
    // TODO logsumexp
    // return raw value for logsumexp
    // return exp(sum);
    return sum;
  }

  double logsumexp(double x, double y, bool init_flag){
    if(init_flag){
      return y;
    }
    if(x == y){
      return x + log(2);
    }
    double vmin = std::min(x, y);
    double vmax = std::max(x, y);
    if(vmax > vmin + 50){
      return vmax;
    }else{
      return vmax + log(exp(vmin - vmax) + 1.0);
    }
  }

  // Calc p(x_t|mu, Sigma) for all topics
  void fit_gaussian(int doc_id){
    // update p(x_t|mu, Sigma)
    std::vector<double> raw_values;
    for(int topic_id = 0; topic_id < this->k_max; ++topic_id){
      raw_values.push_back(this->pi[std::make_pair(topic_id, doc_id - 1)] * calc_prob_doc(doc_id, this->prev_mu[topic_id], this->prev_large_sigma[topic_id]));
    }
    
    double sum_logsumexp = 0.0;
    for(int topic_id = 0; topic_id < this->k_max; ++topic_id){
      sum_logsumexp = logsumexp(sum_logsumexp, raw_values.at(topic_id), (topic_id == 0));
    }

    for(int topic_id = 0; topic_id < this->k_max; ++topic_id){
      if(sum_logsumexp == 0){
	this->prob_doc[std::make_pair(doc_id, topic_id)] = 0.0;
      }else{
	double val = this->pi[std::make_pair(topic_id, doc_id - 1)] * calc_prob_doc(doc_id, this->prev_mu[topic_id], this->prev_large_sigma[topic_id]);
	// this->prob_doc[std::make_pair(doc_id, topic_id)] = log(exp(val)) / sum_logsumexp;
	this->prob_doc[std::make_pair(doc_id, topic_id)] = exp(val - sum_logsumexp);
      }
    }
  }

  void calc_paramaters(int doc_id, int topic_id){
    // update p(i|x_t)
    double t_diff = timestamps[doc_id] - timestamps[doc_id - 1];
    double sum_p = 0.0;
    std::unordered_map<int, double> doc = documents.at(doc_id);
    prob_topic[std::make_pair(topic_id, doc_id)] = pi[std::make_pair(topic_id, doc_id - 1)] * prob_doc[std::make_pair(doc_id, topic_id)];
    
    for(int k = 0; k < k_max; ++k){
      sum_p += pi[std::make_pair(k, doc_id - 1)] * prob_doc[std::make_pair(doc_id, k)];
    }

    if(sum_p == 0){
      prob_topic[std::make_pair(topic_id, doc_id)] = 0.0;
    }else{
      prob_topic[std::make_pair(topic_id, doc_id)] /= sum_p;
    }
    
    // update \gamma_i^(t)
    gamma[std::make_pair(topic_id, doc_id)] = weighted_average(prob_topic[std::make_pair(topic_id, doc_id)],
							       1.0 / this->k_max, 1.0, this->alpha);

    // update \pi_i^(t)
    // in paper, 
    // wrong: \pi_i^(t)
    // correct: \pi_i^(t - 1)
    pi[std::make_pair(topic_id, doc_id)] = weighted_average(pi[std::make_pair(topic_id, doc_id - 1)],
							    gamma[std::make_pair(topic_id, doc_id)],
							    this->prev_m, pow(this->lambda, -(t_diff)));

    // update \mu_i^(t)
    this->now_mu[topic_id] = weighted_average(this->prev_mu[topic_id], doc,
					      pi[std::make_pair(topic_id, doc_id - 1)] * this->prev_m,
					      pow(this->lambda, -t_diff) * gamma[std::make_pair(topic_id, doc_id)]);
    
    // \Gamma_i^t, \Sigma contain only diag elements.
    // update \Gamma_i^t
    this->now_large_gamma[topic_id] = weighted_average(this->prev_large_gamma[topic_id], sqrt(doc),
						       pi[std::make_pair(topic_id, doc_id - 1)] * this->prev_m,
						       pow(this->lambda, -t_diff) * gamma[std::make_pair(topic_id, doc_id)]);

    // update \Sigma
    this->now_large_sigma[topic_id] = sum(this->now_large_gamma[topic_id], prod(-1.0, sqrt(this->now_mu[topic_id])));
  }

  void update_parameter(int doc_id){
    double t_diff = timestamps[doc_id] - timestamps[doc_id - 1];
    this->prev_m = pow(this->lambda, -t_diff) * this->prev_m + 1;
    // burst to inf(pow)
    if(std::isinf(this->prev_m)){
      this->prev_m = 1.0;
    }
    
    for(int topic_id = 0; topic_id < this->k_max; ++topic_id){
      // update prev variable
      this->prev_mu[topic_id] = this->now_mu[topic_id];
      this->prev_large_gamma[topic_id] = this->now_large_gamma[topic_id];
      this->prev_large_sigma[topic_id] = this->now_large_sigma[topic_id];
    }
  }

  void set_initial_parameters(){
    // initialize each params
    // m, t
    this->prev_m = 1.0;
    this->timestamps[-1] = this->timestamps[0];

    for(int topic_id = 0; topic_id < this->k_max; ++topic_id){
      // \pi
      this->pi[std::make_pair(topic_id, -1)] = 1.0 / this->k_max;

      // \Sigma, \mu = rand
      for(int word_id = 0; word_id < words.size(); ++word_id){
	this->prev_large_sigma[topic_id][word_id] = (double)rand() * (1.0 / (RAND_MAX + 1.0));
	this->prev_mu[topic_id][word_id] = (double)rand() * (1.0 / (RAND_MAX + 1.0));
      }
    }

    // scale timestamp to [0, 1]
    double min = this->timestamps[0];
    double max = this->timestamps[documents.size() - 1];
    std::unordered_map<int, double>::iterator i;
    for(i = this->timestamps.begin();i != this->timestamps.end(); ++i){
      (*i).second = ((*i).second - min) / (max - min);
    }
  }

  void fitting(){
    set_initial_parameters();
    for(int doc_id = 0; doc_id < documents.size(); ++doc_id){
      fit_gaussian(doc_id);
      for(int topic_id = 0; topic_id < this->k_max; ++topic_id){
	calc_paramaters(doc_id, topic_id);
      }
      update_parameter(doc_id);
    }
  }
  
private:
  std::vector<std::unordered_map<int, double> > documents;
  std::vector<std::string> words;
  std::unordered_map<int, double> timestamps;

  // key: <doc_id, topic_id>
  /// p(doc_id | topic_id)
  std::unordered_map<key, double, myhash, myeq> prob_doc;
    
  // key: <topic_id, doc_id>
  /// p(topic_id | doc_id)
  std::unordered_map<key, double, myhash, myeq> prob_topic;

  // key: <topic_id, doc_id(= time index t)>
  // value: value
  std::unordered_map<key, double, myhash, myeq> gamma;
  std::unordered_map<key, double, myhash, myeq> pi;

  // update each iterations
  std::unordered_map<int, std::unordered_map<int, double> > now_mu;
  std::unordered_map<int, std::unordered_map<int, double> > now_large_gamma;
  std::unordered_map<int, std::unordered_map<int, double> > now_large_sigma;
  double prev_m;
  std::unordered_map<int, std::unordered_map<int, double> > prev_mu;
  std::unordered_map<int, std::unordered_map<int, double> > prev_large_gamma;
  std::unordered_map<int, std::unordered_map<int, double> > prev_large_sigma;

  // parameters
  double alpha;
  double lambda;
  int k_max;
};

