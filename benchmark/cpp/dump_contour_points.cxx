// Dump full C++ MnContours point coordinates per case
#include "Minuit2/FCNBase.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnStrategy.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/ContoursError.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>
using namespace ROOT::Minuit2;

class Rosen2: public FCNBase { public:
    double operator()(const std::vector<double>& p) const override {
        return (1-p[0])*(1-p[0]) + 100*(p[1]-p[0]*p[0])*(p[1]-p[0]*p[0]);
    } double Up() const override { return 1.0; } };
class Quad4: public FCNBase { public:
    double operator()(const std::vector<double>& p) const override {
        double s=0; for(int i=0;i<4;++i)s+=p[i]*p[i]; return s;
    } double Up() const override { return 1.0; } };
class RosenN: public FCNBase { unsigned N; public:
    RosenN(unsigned n):N(n){}
    double operator()(const std::vector<double>& p) const override {
        double s=0; for(unsigned i=0;i+1<N;++i){double a=p[i],b=p[i+1];
        s+=100*(b-a*a)*(b-a*a)+(1-a)*(1-a);} return s;
    } double Up() const override { return 1.0; } };
class GaussLL: public FCNBase { std::vector<double> d; public:
    GaussLL(){ std::mt19937_64 r(0xCAFEF00D); std::normal_distribution<double> g(2,1);
        d.reserve(100); for(int i=0;i<100;++i)d.push_back(g(r)); }
    double operator()(const std::vector<double>& p) const override {
        double mu=p[0],s=p[1]; if(s<=0) return 1e30; double sum=0;
        for(auto x:d){double dd=x-mu; sum+=std::log(s)+0.5*dd*dd/(s*s);} return sum;
    } double Up() const override { return 0.5; } };
class GaussLLN: public FCNBase { unsigned N; std::vector<double> truths;
    std::vector<std::vector<double>> data; public:
    GaussLLN(unsigned n, unsigned ne):N(n){
        std::mt19937_64 r(0xCAFEF00D); std::normal_distribution<double> g(0,1);
        truths.resize(n); for(unsigned i=0;i<n;++i) truths[i]=g(r);
        unsigned per=std::max(1u, ne/n);
        data.assign(n, std::vector<double>(per));
        for(unsigned i=0;i<n;++i) for(unsigned j=0;j<per;++j) data[i][j]=truths[i]+g(r);
    }
    double operator()(const std::vector<double>& p) const override {
        double s=0; for(unsigned i=0;i<N;++i){double mu=p[i]; for(double x:data[i])
        {double d=x-mu; s+=0.5*d*d;}} return s;
    } double Up() const override { return 0.5; } };

template<class F> void dump(const std::string& name, F&& fcn,
                              const std::vector<std::pair<std::string,double>>& init) {
    MnUserParameters upar;
    for (auto& [n, v] : init) upar.Add(n.c_str(), v, 0.1);
    MnMigrad migrad(fcn, upar, MnStrategy(0));
    FunctionMinimum mn = migrad();
    MnContours contours(fcn, mn, MnStrategy(0));
    ContoursError ce = contours.Contour(0, 1, 30);
    std::cout << "== " << name << " ==\n";
    std::cout << "  npts=" << ce().size() << "  nfcn=" << ce.NFcn() << "\n";
    std::cout << "  fmin_x=" << std::setprecision(15) << mn.UserState().Value(0)
              << "  fmin_y=" << mn.UserState().Value(1) << "\n";
    for (size_t i = 0; i < ce().size(); ++i) {
        std::cout << "  pt[" << i << "] = (" << std::setprecision(15)
                  << ce()[i].first << ", " << ce()[i].second << ")\n";
    }
}

int main() {
    dump("rosenbrock_2d",  Rosen2{}, {{"p0",-1.2},{"p1",1.0}});
    dump("rosenbrock_10d", RosenN(10),
         {{"p0",-1.2},{"p1",1.0},{"p2",-1.2},{"p3",1.0},{"p4",-1.2},
          {"p5",1.0},{"p6",-1.2},{"p7",1.0},{"p8",-1.2},{"p9",1.0}});
    dump("quad_4d",        Quad4{}, {{"p0",1.0},{"p1",1.0},{"p2",1.0},{"p3",1.0}});
    dump("gauss_ll_2_100", GaussLL{}, {{"mu",1.0},{"sigma",2.0}});
    dump("gauss_ll_10_1000", GaussLLN(10,1000),
         {{"p0",0.0},{"p1",0.0},{"p2",0.0},{"p3",0.0},{"p4",0.0},
          {"p5",0.0},{"p6",0.0},{"p7",0.0},{"p8",0.0},{"p9",0.0}});
    return 0;
}
