#include<iostream>
#include <fstream>
#include<cmath>
#include<math.h>
#include<vector>
#include <iomanip>
#include<string>
int iter = 0;
int culc_func = 0;
const double alpha1 = 2.0;
const double alpha2 = 17.0;
int Constobn =2;
int indf = 3;
int h1 = 0;
int h2 = 0;
int count = 0;
double p, l;
double AimFunc(const double x, const double y)
{
    culc_func++;
    if (indf == 1)
    {
        return 10 * x * x - 4 * x * y + 7 * y * y - 4 * sqrt(5) * (5 * x - y) - 16;
    }
    if (indf == 2)
    {
        return alpha1 * (x * x - y) * (x * x - y) + (x - 1) * (x - 1);
    }
    if (indf == 3)
    {
        return alpha2 * (x * x - y) * (x * x - y) + (x - 1) * (x - 1);
    }
}
double Funcx(const double x, const double y)
{
    culc_func++;
    if (indf == 1)
    {
        return 20 * x - 4 * y - 20 * sqrt(5);
    }
    if (indf == 2)
    {
        return -2 + 2 * x + 4 * x * x * x * alpha1 - 4 * x * y * alpha1;
    }
    if (indf == 3)
    {
        return -2 + 2 * x + 4 * x * x * x * alpha2 - 4 * x * y * alpha2;
    }
}
double Funcy(const double x, const double y)
{
    culc_func++;
    if (indf == 1)
    {
        return  4 * sqrt(5) - 4 * x + 14 * y;
    }
    if (indf == 2)
    {
        return -2 * x * x * alpha1 + 2 * y * alpha1;
    }
    if (indf == 3)
    {
        return -2 * x * x * alpha2 + 2 * y * alpha2;
    }
}
double AimFunclymbda(const double lymbda, const double pointx, const double pointy, const double fx, const double fy)
{
    double x, y;
    x = pointx + lymbda * fx;
    y = pointy + lymbda * fy;
    return AimFunc(x, y);
}
//double Funcxx(const double x, const double y)
//{
//    culc_func++;
//    if (indf == 1)
//    {
//        return 20;
//    }
//    if (indf == 2)
//    {
//        return 2 + 12 * x * x * alpha1 - 4 * y * alpha1;
//    }
//    if (indf == 3)
//    {
//        return 2 + 12 * x * x * alpha2 - 4 * y * alpha2;
//    }
//}
//double Funcyy(const double x, const double y)
//{
//    culc_func++;
//    if (indf == 1)
//    {
//        return 14;
//    }
//    if (indf == 2)
//    {
//        return 2 * alpha1;
//    }
//    if (indf == 3)
//    {
//        return 2 * alpha2;
//    }
//}
//double Funcxy(const double x, const double y)
//{
//    culc_func++;
//    if (indf == 1)
//    {
//        return -4;
//    }
//    if (indf == 2)
//    {
//        return -4 * x * alpha1;
//    }
//    if (indf == 3)
//    {
//        return -4 * x * alpha2;
//    }
//}
//double Pathx(double  fxx, double fxy, double fyy, double fx, double fy)
//{
//    return (fx * fyy - fy * fxy) / (fxx * fyy - fxy * fxy);
//}
//double Pathy(double  fxx, double fxy, double fyy, double fx, double fy)
//{
//    return (fxx * fy - fxy * fx) / (fxx * fyy - fxy * fxy);
//}
//double NormVecGramm(double x, double y, double xt, double yt)
//{
//    return sqrt(NormVec(Funcx(xt, yt), Funcx(xt, yt)) * x * x + y * y);
//}
double NormVec(double x, double y)
{
    return sqrt(x * x + y * y);
}
double Deltawdeltap(double deltawx, double deltawy, double deltapx, double deltapy)
{
    return (deltawx * deltapx + deltawy * deltapy);
}
double DeltawA(double a11, double a12, double a21, double a22, double deltawx, double deltawy)
{
    return (a11 * deltawx + a12 * deltawy) * deltawx + (a21 * deltawx + a22 * deltawy) * deltawy;
}
double PathAMGx(double a11, double a12, double fx, double fy)
{
    return a11 * fx + a12 * fy;
}
double PathAMGy(double a21, double a22, double fx, double fy)
{
    return a21 * fx + a22 * fy;
}
double deltaA11DFP(double a11, double a12, double a21, double a22, double deltawx, double deltawy, double deltapx, double deltapy)
{
    double a = (deltapx * deltapx) / Deltawdeltap(deltawx, deltawy, deltapx, deltapy);
    double b = (a11*a11*deltawx*deltawx+a12*a12* deltawy * deltawy+2*a11*a12* deltawx * deltawy) / DeltawA(a11, a12, a21, a22, deltawx, deltawy);
    return -(a + b);
}
double deltaA12DFP(double a11, double a12, double a21, double a22, double deltawx, double deltawy, double deltapx, double deltapy)
{
    double a = (deltapx * deltapy) / Deltawdeltap(deltawx, deltawy, deltapx, deltapy);
    double b = (a11*a21* deltawx * deltawx + a12*a22* deltawy * deltawy + a12*a21* deltawx * deltawy + a11*a22* deltawx * deltawy) / DeltawA(a11, a12, a21, a22, deltawx, deltawy);
    return -(a + b);
}
double deltaA21DFP(double a11, double a12, double a21, double a22, double deltawx, double deltawy, double deltapx, double deltapy)
{
    double a = (deltapx * deltapy) / Deltawdeltap(deltawx, deltawy, deltapx, deltapy);
    double b = (a11 * a21 * deltawx * deltawx + a12 * a22 * deltawy * deltawy + a12 * a21 * deltawx * deltawy + a11 * a22 * deltawx * deltawy) / DeltawA(a11, a12, a21, a22, deltawx, deltawy);
    return -(a + b);
}
double deltaA22DFP(double a11, double a12, double a21, double a22, double deltawx, double deltawy, double deltapx, double deltapy)
{
    double a = (deltapy * deltapy) / Deltawdeltap(deltawx, deltawy, deltapx, deltapy);
    double b = (a21*a21*deltawx*deltawx+a22*a22*deltawy * deltawy+2*a21*a22*deltawx * deltawy) / DeltawA(a11, a12, a21, a22, deltawx, deltawy);
    return -(a + b);
}
double deltaA11BFSH(double a11, double a12, double a21, double a22, double deltawx, double deltawy, double deltapx, double deltapy)
{
    double a = (deltapx * deltapx) / Deltawdeltap(deltawx, deltawy, deltapx, deltapy);
    double b = (a11 * a11 * deltawx * deltawx + a12 * a12 * deltawy * deltawy + 2 * a11 * a12 * deltawx * deltawy) / DeltawA(a11, a12, a21, a22, deltawx, deltawy);
    double c = DeltawA(a11, a12, a21, a22, deltawx, deltawy) *((a11 * deltawx + a12 * deltawy) / DeltawA(a11, a12, a21, a22, deltawx, deltawy) - deltapx / Deltawdeltap(deltawx, deltawy, deltapx, deltapy))* ((a11 * deltawx + a12 * deltawy) / DeltawA(a11, a12, a21, a22, deltawx, deltawy) - deltapx / Deltawdeltap(deltawx, deltawy, deltapx, deltapy));
    return -(a + b)+c;
}
double deltaA12BFSH(double a11, double a12, double a21, double a22, double deltawx, double deltawy, double deltapx, double deltapy)
{
    double a = (deltapx * deltapy) / Deltawdeltap(deltawx, deltawy, deltapx, deltapy);
    double b = (a11 * a21 * deltawx * deltawx + a12 * a22 * deltawy * deltawy + a12 * a21 * deltawx * deltawy + a11 * a22 * deltawx * deltawy) / DeltawA(a11, a12, a21, a22, deltawx, deltawy);
    double c = DeltawA(a11, a12, a21, a22, deltawx, deltawy) * ((a11 * deltawx + a12 * deltawy) / DeltawA(a11, a12, a21, a22, deltawx, deltawy) - deltapx / Deltawdeltap(deltawx, deltawy, deltapx, deltapy)) * ((a21 * deltawx + a22 * deltawy) / DeltawA(a11, a12, a21, a22, deltawx, deltawy) - deltapy / Deltawdeltap(deltawx, deltawy, deltapx, deltapy));
    return -(a + b)+c;
}
double deltaA21BFSH(double a11, double a12, double a21, double a22, double deltawx, double deltawy, double deltapx, double deltapy)
{
    double a = (deltapx * deltapy) / Deltawdeltap(deltawx, deltawy, deltapx, deltapy);
    double b = (a11 * a21 * deltawx * deltawx + a12 * a22 * deltawy * deltawy + a12 * a21 * deltawx * deltawy + a11 * a22 * deltawx * deltawy) / DeltawA(a11, a12, a21, a22, deltawx, deltawy);
    double c = DeltawA(a11, a12, a21, a22, deltawx, deltawy) * ((a11 * deltawx + a12 * deltawy) / DeltawA(a11, a12, a21, a22, deltawx, deltawy) - deltapx / Deltawdeltap(deltawx, deltawy, deltapx, deltapy)) * ((a21 * deltawx + a22 * deltawy) / DeltawA(a11, a12, a21, a22, deltawx, deltawy) - deltapy / Deltawdeltap(deltawx, deltawy, deltapx, deltapy));
    return -(a + b)+c;
}
double deltaA22BFSH(double a11, double a12, double a21, double a22, double deltawx, double deltawy, double deltapx, double deltapy)
{
    double a = (deltapy * deltapy) / Deltawdeltap(deltawx, deltawy, deltapx, deltapy);
    double b = (a21 * a21 * deltawx * deltawx + a22 * a22 * deltawy * deltawy + 2 * a21 * a22 * deltawx * deltawy) / DeltawA(a11, a12, a21, a22, deltawx, deltawy);
    double c = DeltawA(a11, a12, a21, a22, deltawx, deltawy) * ((a21 * deltawx + a22 * deltawy) / DeltawA(a11, a12, a21, a22, deltawx, deltawy) - deltapy / Deltawdeltap(deltawx, deltawy, deltapx, deltapy)) * ((a21 * deltawx + a22 * deltawy) / DeltawA(a11, a12, a21, a22, deltawx, deltawy) - deltapy / Deltawdeltap(deltawx, deltawy, deltapx, deltapy));
    return -(a + b)+c;
}
double deltaA11P(double a11, double a12, double a21, double a22, double deltawx, double deltawy, double deltapx, double deltapy)
{
    return (deltapx-a11*deltawx-a12*deltawy)* (deltapx - a11 * deltawx - a12 * deltawy)/((deltapx - a11 * deltawx - a12 * deltawy)*deltawx+(deltapy-a21*deltawx-a22*deltawy)*deltawy);
}
double deltaA12P(double a11, double a12, double a21, double a22, double deltawx, double deltawy, double deltapx, double deltapy)
{
    return (deltapx - a11 * deltawx - a12 * deltawy) * (deltapy - a21 * deltawx - a22 * deltawy) / ((deltapx - a11 * deltawx - a12 * deltawy) * deltawx + (deltapy - a21 * deltawx - a22 * deltawy) * deltawy);
}
double deltaA21P(double a11, double a12, double a21, double a22, double deltawx, double deltawy, double deltapx, double deltapy)
{
    return (deltapx - a11 * deltawx - a12 * deltawy) * (deltapy - a21 * deltawx - a22 * deltawy) / ((deltapx - a11 * deltawx - a12 * deltawy) * deltawx + (deltapy - a21 * deltawx - a22 * deltawy) * deltawy);
}
double deltaA22P(double a11, double a12, double a21, double a22, double deltawx, double deltawy, double deltapx, double deltapy)
{
    return (deltapy - a21 * deltawx - a22 * deltawy) * (deltapy - a21 * deltawx - a22 * deltawy) / ((deltapx - a11 * deltawx - a12 * deltawy) * deltawx + (deltapy - a21 * deltawx - a22 * deltawy) * deltawy);
}
double Gold(const double eps, const double pointx, const double pointy, const double fx, const double fy)
{
    const double tau1 = 2 / (sqrt(5) + 3);
    const double tau2 = 2 / (sqrt(5) + 1);
    double x1, x2, f1, f2;
    double ak = 0;
    double bk;
    if (iter >=3) 
    {
        bk = l / NormVec(fx, fy);
    }
    else
    {
        bk = p / NormVec(fx, fy);
    }
    x1 = ak + tau1 * (bk - ak);
    x2 = ak + tau2 * (bk - ak);
    f1 = AimFunclymbda(x1, pointx, pointy, fx, fy);
    f2 = AimFunclymbda(x2, pointx, pointy, fx, fy);
    while (bk - ak > eps*eps)
    {
        if (f1 < f2)
        {
            bk = x2;
            x2 = x1;
            f2 = f1;
            x1 = ak + tau1 * (bk - ak);
            f1 = AimFunclymbda(x1, pointx, pointy, fx, fy);
        }
        else
        {
            ak = x1;
            x1 = x2;
            f1 = f2;
            x2 = ak + tau2 * (bk - ak);
            f2 = AimFunclymbda(x2, pointx, pointy, fx, fy);
        }
    }
    return (bk+ak)/2;
}
std::vector<double> DFP(double eps, double xn, double yn)
{
    std::vector<double> points;
    double lymbda;
    double deltaa11,deltaa12,deltaa21,deltaa22,a11,a12,a21,a22,xt, yt,xtl,ytl,fxl,fyl, fx, fy, px, py,deltax,deltay,deltawx,deltawy;
    int indh = 0;
    xt = xn;
    yt = yn;
    a11 = 1.0;
    a12 = 0.0;
    a21 = 0.0;
    a22 = 1.0;
    fx = -Funcx(xt, yt);
    fy = -Funcy(xt, yt);
    points.push_back(xt);
    points.push_back(yt);
    while (NormVec(fx, fy) > eps)
    {
        iter++;
        count++;
        px = PathAMGx(a11, a21, fx, fy);
        py = PathAMGy(a12, a22, fx, fy);
        lymbda = Gold(eps, xt, yt, px, py);
        xtl = xt;
        ytl = yt;
        xt = xt + lymbda * px;
        yt = yt + lymbda * py;
        deltax = xt-xtl;
        deltay = yt-ytl;
        fxl = fx;
        fyl = fy;
        fx = -Funcx(xt, yt);
        fy = -Funcy(xt, yt);
        deltawx = fx- fxl;
        deltawy = fy - fyl;
        deltaa11= deltaA11DFP(a11, a12, a21, a22, deltawx, deltawy, deltax, deltay);
        deltaa12= deltaA12DFP(a11, a12, a21, a22, deltawx, deltawy, deltax, deltay);
        deltaa21= deltaA21DFP(a11, a12, a21, a22, deltawx, deltawy, deltax, deltay);
        deltaa22= deltaA22DFP(a11, a12, a21, a22, deltawx, deltawy, deltax, deltay);
        a11 = a11 + deltaa11;
        a12 = a12 + deltaa12;
        a21 = a21 + deltaa21;
        a22 = a22 + deltaa22;
        points.push_back(xt);
        points.push_back(yt);
        if (count == Constobn)
        {
            count = 0;
            a11 = 1.;
            a12 = 0.;
            a21 = 0.;
            a22 = 1.;
        }
    }
    return points;
}
std::vector<double> BFSH(double eps, double xn, double yn)
{
    std::vector<double> points;
    double lymbda;
    double deltaa11, deltaa12, deltaa21, deltaa22, a11, a12, a21, a22, xt, yt, xtl, ytl, fxl, fyl, fx, fy, px, py, deltax, deltay, deltawx, deltawy;
    int indh = 0;
    xt = xn;
    yt = yn;
    a11 = 1.0;
    a12 = 0.0;
    a21 = 0.0;
    a22 = 1.0;
    fx = -Funcx(xt, yt);
    fy = -Funcy(xt, yt);
    points.push_back(xt);
    points.push_back(yt);
    while (NormVec(fx, fy) > eps)
    {
        iter++;
        count++;
        px = PathAMGx(a11, a12, fx, fy);
        py = PathAMGy(a21, a22, fx, fy);
        lymbda = Gold(eps, xt, yt, px, py);
        xtl = xt;
        ytl = yt;
        xt = xt + lymbda * px;
        yt = yt + lymbda * py;
        deltax = xt - xtl;
        deltay = yt - ytl;
        fxl = fx;
        fyl = fy;
        fx = -Funcx(xt, yt);
        fy = -Funcy(xt, yt);
        deltawx = fx - fxl;
        deltawy = fy - fyl;
        deltaa11 = deltaA11BFSH(a11, a12, a21, a22, deltawx, deltawy, deltax, deltay);
        deltaa12 = deltaA12BFSH(a11, a12, a21, a22, deltawx, deltawy, deltax, deltay);
        deltaa21 = deltaA21BFSH(a11, a12, a21, a22, deltawx, deltawy, deltax, deltay);
        deltaa22 = deltaA22BFSH(a11, a12, a21, a22, deltawx, deltawy, deltax, deltay);
        a11 = a11 + deltaa11;
        a12 = a12 + deltaa12;
        a21 = a21 + deltaa21;
        a22 = a22 + deltaa22;
        points.push_back(xt);
        points.push_back(yt);
        if (count == Constobn)
        {
            count = 0;
            a11 = 1.;
            a12 = 0.;
            a21 = 0.;
            a22 = 1.;
        }
        /*std::cout << a11 << ' ' << a12 << std::endl;
        std::cout << a21 << ' ' << a22 << std::endl;*/
    }
    return points;
}
std::vector<double> PAUL(double eps, double xn, double yn)
{
    std::vector<double> points;
    double lymbda = 1;
    double deltaa11, deltaa12, deltaa21, deltaa22, a11, a12, a21, a22, xt, yt, xtl, ytl, fxl, fyl, fx, fy, px, py, deltax, deltay, deltawx, deltawy;
    int indh = 0;
    xt = xn;
    yt = yn;
    a11 = 1.0;
    a12 = 0.0;
    a21 = 0.0;
    a22 = 1.0;
    fx = -Funcx(xt, yt);
    fy = -Funcy(xt, yt);
    points.push_back(xt);
    points.push_back(yt);
    while (NormVec(fx, fy) > eps)
    {
        iter++;
        count++;
        px = PathAMGx(a11, a12, fx, fy);
        py = PathAMGy(a21, a22, fx, fy);
        lymbda = Gold(eps, xt, yt, px, py);
        if (lymbda < 10*eps)
        {
            count = 0;
            a11 = 1.;
            a12 = 0.;
            a21 = 0.;
            a22 = 1.;
            px = PathAMGx(a11, a12, fx, fy);
            py = PathAMGy(a21, a22, fx, fy);
            lymbda = Gold(eps, xt, yt, px, py);
        }
        xtl = xt;
        ytl = yt;
        xt = xt + lymbda * px;
        yt = yt + lymbda * py;
        deltax = xt - xtl;
        deltay = yt - ytl;
        fxl = fx;
        fyl = fy;
        fx = -Funcx(xt, yt);
        fy = -Funcy(xt, yt);
        deltawx = fx - fxl;
        deltawy = fy - fyl;
        deltaa11 = deltaA11P(a11, a12, a21, a22, deltawx, deltawy, deltax, deltay);
        deltaa12 = deltaA12P(a11, a12, a21, a22, deltawx, deltawy, deltax, deltay);
        deltaa21 = deltaA21P(a11, a12, a21, a22, deltawx, deltawy, deltax, deltay);
        deltaa22 = deltaA22P(a11, a12, a21, a22, deltawx, deltawy, deltax, deltay);
        a11 = a11 + deltaa11;
        a12 = a12 + deltaa12;
        a21 = a21 + deltaa21;
        a22 = a22 + deltaa22;
        points.push_back(xt);
        points.push_back(yt);
        if (lymbda<eps)
        {
            count = 0;
            a11 = 1.;
            a12 = 0.;
            a21 = 0.;
            a22 = 1.;
        }
        /*std::cout << a11 << ' ' << a12 << std::endl;
        std::cout << a21 << ' ' << a22 << std::endl;
        std::cout <<  std::endl;
        std::cout << xt << ' ' << yt << std::endl;
        std::cout << std::endl;
        std::cout << fx << ' ' << fy<< std::endl;
        std::cout << std::endl;
        std::cout << px << ' ' << py <<' '<< lymbda<< std::endl;
        std::cout << std::endl;*/
    }
    return points;
}
//std::vector<double> NewtonStep(double eps, double xn, double yn)
//{
//    std::vector<double> points;
//    double lymbda = 0;
//    double xt, yt, fx, fy, fxx, fxy, fyy, px, py;
//    double omegat;
//    double etta = 0;
//    double dzetta = 0;
//    int indh = 0;
//    xt = xn;
//    yt = yn;
//    fx = -Funcx(xt, yt);
//    fy = -Funcy(xt, yt);
//    points.push_back(xt);
//    points.push_back(yt);
//    while (NormVec(fx, fy) > eps)
//    {
//        indh = 0;
//        iter++;
//        fxx = Funcxx(xt, yt);
//        fxy = Funcxy(xt, yt);
//        fyy = Funcyy(xt, yt);
//        if ((fxx > 0) && ((fxx * fyy - fxy * fxy) > 0))
//        {
//            indh = 0;
//            px = Pathx(fxx, fxy, fyy, fx, fy);
//            py = Pathy(fxx, fxy, fyy, fx, fy);
//        }
//        else
//        {
//            indh = 1;
//            while (indh == 1)
//            {
//                fxx = fxx + eps;
//                fyy = fyy + eps;
//                if ((fxx > 0) && ((fxx * fyy - fxy * fxy) > 0))
//                {
//                    h2++;
//                    indh = 0;
//                }
//            }
//            px = Pathx(fxx, fxy, fyy, fx, fy);
//            py = Pathy(fxx, fxy, fyy, fx, fy);
//        }
//        lymbda = Gold(eps, xt, yt, px, py);
//        xt = xt + lymbda * px;
//        yt = yt + lymbda * py;
//        fx = -Funcx(xt, yt);
//        fy = -Funcy(xt, yt);
//        points.push_back(xt);
//        points.push_back(yt);
//    }
//    points.push_back(xt);
//    points.push_back(yt);
//
//    return points;
//}

int main()
{
    const double eps[2] = { 0.01,0.00001 };
    std::vector<double> points1, points2;
    std::ofstream file;
    std::string fil, end;
    double x0[3] = { 2,1.5,-1 };
    double y0[3] = { -5,2,10 };
    p =10;
    l = 2;
    Constobn = 100;  // Константа обновления
    indf = 2;   //функция
    int i = 1;   //точность
    int j = 3; //точка
    fil = "DFPeps";
    end = std::to_string(i);
    end += "f";
    end += std::to_string(indf);
    end += "p";
    end += std::to_string(j);
    end += ".txt";
    std::cout << std::setprecision(7);
    file << std::setprecision(7);
    file.open(fil + end);
    points1 = DFP(eps[i - 1], x0[j - 1], y0[j - 1]);
    file << iter << ' ' << culc_func << std::endl;
    file << std::setprecision(7);
    for (int i = 0; i < points1.size() - 1; i += 2)
    {
        file << points1[i] << ' ' << points1[i + 1] << std::endl;
    }
    file.close();
    std::cout << "DFP:" << std::endl;
    std::cout << iter << ' ' << culc_func << std::endl;
    std::cout << std::endl;
    culc_func = 0;
    iter = 0;
    count = 0;
    fil = "BFSHeps";
    file.open(fil + end);
    points2 = BFSH(eps[i - 1], x0[j - 1], y0[j - 1]);
    file << iter << ' ' << culc_func << std::endl;
    for (int i = 0; i < points2.size(); i += 2)
    {
        file << points2[i] << ' ' << points2[i + 1] << std::endl;
    }
    file.close();
    std::cout << "BFSH:" << std::endl;
    std::cout << iter << ' ' << culc_func << std::endl;
    std::cout << std::endl;
    culc_func = 0;
    iter = 0;
    count = 0;
    fil = "PAULeps";
    file.open(fil + end);
    points2 = PAUL(eps[i - 1], x0[j - 1], y0[j - 1]);
    file << iter << ' ' << culc_func << std::endl;
    for (int i = 0; i < points2.size(); i += 2)
    {
        file << points2[i] << ' ' << points2[i + 1] << std::endl;
    }
    file.close();
    std::cout << "PAUL:" << std::endl;
    std::cout << iter << ' ' << culc_func << std::endl;
    std::cout << std::endl;
    std::cout << "Hello World!\n";
    return 0;
}


