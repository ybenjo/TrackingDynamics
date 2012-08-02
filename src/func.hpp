#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <unordered_map>
#include <vector>
#include <string>
#include <algorithm>
#include <stdlib.h>

#define _USE_MATH_DEFINES
#include <cmath>

struct myeq : std::binary_function<std::pair<int, int>, std::pair<int, int>, bool>{
  bool operator() (const std::pair<int, int> & x, const std::pair<int, int> & y) const{
    return x.first == y.first && x.second == y.second;
  }
};

struct myhash : std::unary_function<std::pair<int, int>, size_t>{
private:
  const std::hash<int> h_int;
public:
  myhash() : h_int() {}
  size_t operator()(const std::pair<int, int> & p) const{
    size_t seed = h_int(p.first);
    return h_int(p.second) + 0x9e3779b9 + (seed<<6) + (seed>>2);
  }
};

std::vector<std::string> split_string(std::string s, std::string c){
  std::vector<std::string> ret;
  for(int i = 0, n = 0; i <= s.length(); i = n + 1){
    n = s.find_first_of(c, i);
    if(n == std::string::npos) n = s.length();
    std::string tmp = s.substr(i, n-i);
    ret.push_back(tmp);
  }
  return ret;
}

std::unordered_map<int, double> prod(double x, const std::unordered_map<int, double>& y){
  std::unordered_map<int, double> ret;
  std::unordered_map<int, double>::const_iterator i;
  for(i = y.begin(); i != y.end(); ++i){
    double val = (*i).second;
    ret[(*i).first] = x * val;
  }
  return ret;
}

std::unordered_map<int, double> sum(const std::unordered_map<int, double>& x, const std::unordered_map<int, double>& y){
  std::unordered_map<int, double> ret;
  std::unordered_map<int, double>::const_iterator i;
  for(i = x.begin(); i != x.end(); ++i){
    int key = (*i).first;
    ret[key] += (*i).second;
  }
  
  for(i = y.begin(); i != y.end(); ++i){
    int key = (*i).first;
    ret[key] += (*i).second;
  }
  
  return ret;
}

std::unordered_map<int, double> sqrt(const std::unordered_map<int, double>& x){
  std::unordered_map<int, double> ret;
  std::unordered_map<int, double>::const_iterator i;
  for(i = x.begin(); i != x.end(); ++i){
    ret[(*i).first] = pow((*i).second, 2);
  }
  return ret;
}

double weighted_average(double x, double y, double a, double b){
  return (1 - (b/ (a + b))) * x + (1 - (a/ (a + b))) * y;
}

std::unordered_map<int, double> weighted_average(const std::unordered_map<int, double>& x,
						 const std::unordered_map<int, double>& y,
						 double a, double b){
  std::unordered_map<int, double> ret;
  ret = sum(prod((a / (a + b)), x), prod((b / (a + b)), y));
  return ret;
}

std::unordered_map<int, double> inverse(const std::unordered_map<int, double>& x){
  std::unordered_map<int, double> ret;
  std::unordered_map<int, double>::const_iterator i;
  for(i = x.begin(); i != x.end(); ++i){
    ret[(*i).first] = 1 / (*i).second;
  }
  return ret;
}
