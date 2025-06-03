#include <iostream>
#include <fstream>
#include<cmath>
#include<math.h>
#include<vector>
#include <iomanip>
int iter = 0;
int culc_func = 0;
const double alpha1 = 1.0;
const double alpha2 = 5.0;
double p = 3.5754019;
double AimFunc(const double x,const double y)
{
    culc_func++;
    return alpha2*(x*x-y)*(x*x-y)+(x-1)*(x-1);
   /* return 10*x*x-4*x*y+7*y*y-4*sqrt(5)*(5*x-y)-16;*/
}
double Funcx(const double x,const double y)
{
    culc_func++;
    /*return 20  * x - 4 * y - 20 * sqrt(5);*/
    return 4*alpha2*x*(x*x-y)+2*(x-1);
}
double Funcy(const double x, const double y)
{
    culc_func++;
    /*return - 4 * x  + 14 * y  + 4 * sqrt(5);*/
    return -2*alpha2*(x*x-y);
}
double AimFunclymbda(const double lymbda,const double pointx,const double pointy,const double fx, const double fy)
{
    culc_func++;
    double x, y;
    x = pointx - lymbda * fx;
    y = pointy - lymbda * fy;
    /*return 10 * x * x - 4 * x * y + 7 * y * y - 4 * sqrt(5) * (5 * x - y) - 16;*/
    return alpha2 * (x * x - y) * (x * x - y) + (x - 1) * (x - 1);
}
double NormVec(double x, double y)
{
    return sqrt(x * x + y * y);
}
double Gold(const double eps,const double pointx,const double pointy)
{
    const double tau1 = 2 / (sqrt(5) + 3);
    const double tau2 = 2 / (sqrt(5) + 1);
    double x1,x2,f1,f2,fx,fy;
    fx = Funcx(pointx, pointy);
    fy= Funcy(pointx, pointy);
    double ak = 0;
    if (iter >= 2)
    {
        p = 1;
    }
    double bk = p*(1/NormVec(fx,fy));
   
    int count = 0;
    x1 = ak + tau1 * (bk - ak);
    x2 = ak + tau2 * (bk - ak);
    f1 = AimFunclymbda(x1,pointx,pointy,fx,fy);
    f2 = AimFunclymbda(x2, pointx, pointy,fx,fy);
    while (bk - ak > eps)
    {
        count++;
        if (f1 <= f2)
        {
            bk = x2;
            x2 = x1;
            f2 = f1;
            x1 = ak + tau1 * (bk - ak);
            f1 = AimFunclymbda(x1, pointx, pointy,fx,fy);
        }
        else
        {
            ak = x1;
            x1 = x2;
            f1 = f2;
            x2 = ak + tau2 * (bk - ak);
            f2 = AimFunclymbda(x2, pointx, pointy,fx,fy);
        }
    }
    return (bk + ak) / 2;
}
std::vector<double> Steepest(double eps, double xn, double yn)
{
    std::vector<double> points;
    double lymbda=0;
    double xt, yt,fx,fy;
    double omegat;
    xt = xn;
    yt = yn;
    fx = Funcx(xt, yt);
    fy = Funcy(xt, yt);
    omegat = NormVec(fx, fy);
    while ( omegat> eps/10)
    {
        iter++;
        lymbda = Gold(eps, xt, yt);
        points.push_back(xt);
        points.push_back(yt);
        xt = xt - lymbda * fx;
        yt = yt - lymbda * fy;
        fx = Funcx(xt, yt);
        fy = Funcy(xt, yt);
        omegat = NormVec(fx, fy);
    }
    points.push_back(xt);
    points.push_back(yt);
    
    return points;
}
std::vector<double> Spleeting(double eps,double xn,double yn)
{
    std::vector<double> points;
    double nu;
    double xt, yt,xt1,yt1;
    double kapa=0.1;
    double f1, f2,fx,fy;
    double omegat;
    double omega = 0.5;
    xt = xn;
    yt = yn;
    fx = Funcx(xt, yt);
    fy = Funcy(xt, yt);
    omegat = NormVec(fx,fy);
    while (omegat > eps/10)
    {
        iter++;
        points.push_back(xt);
        points.push_back(yt);
        nu = 0.5;
        kapa = 0.1;
        omegat = NormVec(fx, fy);
        xt1 = xt - kapa * fx;
        yt1 = yt - kapa * fy;
        f1 = AimFunc(xt, yt);
        f2 = AimFunc(xt1, yt1);

        while (f1 - f2 < omega * kapa * omegat * omegat)
        {
            kapa = nu * kapa;
            xt1 = xt - kapa * fx;
            yt1 = yt - kapa * fy;
            f2 = AimFunc(xt1, yt1);
        }
        xt = xt - kapa * fx;
        yt = yt - kapa * fy;
        fx = Funcx(xt, yt);
        fy = Funcy(xt, yt);
        omegat = NormVec(fx, fy);

    }
    points.push_back(xt);
    points.push_back(yt);
    return points;

}
int main()
{
    const double eps1 = 0.01;
    const double eps2 = 0.00001;
    std::vector<double> points1,points2;
    std::ofstream file;
    points1=Spleeting(eps1, 2,2);
    file.open("Spleetingeps1f3p1.txt");
    file << std::setprecision(7);
    for (int i = 0; i < points1.size()-1; i+=2)
    {
        file << points1[i] << ' ' << points1[i + 1] << std::endl;
    }
    file.close();
    std::cout << std::endl;
    std::cout <<iter<<' ' << culc_func << std::endl;
    std::cout<<std::endl;
    culc_func = 0;
    iter = 0;
    file.open("Steepesteps1f3p1.txt");
    points2 = Steepest(eps1, 2,2);
    for (int i = 0; i < points2.size(); i+=2)
    {
        file << points2[i] << ' ' << points2[i + 1] << std::endl;
    }
    file.close();
    std::cout << std::endl;
    std::cout << iter << ' ' << culc_func << std::endl;
    std::cout << std::endl;
    std::cout << "Hello World!\n";
    return 0;
}

