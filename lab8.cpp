#include<iostream>
#include <fstream>
#include<cmath>
#include<math.h>
#include<vector>
#include <iomanip>
#include<string>
int iter = 0;
int iter1 = 0;
int culc_func = 0;
int indf;
int indogr;
double r, betta;
double p, l;
int indu = 0;
int indc = 0;
double Ogr2(const double x, const double y)
{
    if (indogr == 1)
    {
        return  r * (((-x - 10) +fabs(-x-10))/2+ (-y+fabs(-y))/2 + ((-1 + x + y) + fabs(-1 + x + y)) / 2)* (((-x - 10) + fabs(-x - 10)) / 2 + (-y + fabs(-y)) / 2 + ((-1 + x + y) + fabs(-1 + x + y)) / 2);
    }
    else
    {
        return  r * (((x + 3) * (x + 3) / 16 + (y + 4) * (y + 4) / 9 - 1 + fabs((x + 3) * (x + 3) / 16 + (y + 4) * (y + 4) / 9 - 1)) / 2)* (((x + 3) * (x + 3) / 16 + (y + 4) * (y + 4) / 9 - 1 + fabs((x + 3) * (x + 3) / 16 + (y + 4) * (y + 4) / 9 - 1)) / 2);
    }
}
double Ogr(const double x, const double y)
{
    if (indogr == 1)
    {
        return  r * (1 / (x + 10) + 1 / (y)+1 / (1 - x - y));
    }       
    else
    {
        return  r * (1 / (1 - (x + 3) * (x + 3) / 16 - (y + 4) * (y + 4) / 9));
    }
}
double Aim(const double x, const double y)
{
    if (indf == 1)
    {
        return 10 * x * x - 4 * x * y + 7 * y * y - 4 * sqrt(5) * (5 * x - y) - 16;
    }
    if (indf == 2)
    {
        return 1 * (x * x - y) * (x * x - y) + (x - 1) * (x - 1);
    }
    if (indf == 3)
    {
        return 12 * (x * x - y) * (x * x - y) + (x - 1) * (x - 1);
    }
}
double AimFunc(const double x, const double y)
{
    culc_func++;
    return Aim(x, y) + Ogr(x, y);
}
double AimFunc2(const double x, const double y)
{
    culc_func++;
    return Aim(x, y) + Ogr2(x, y);
}
int CritOgr(double x,double y)
{
    culc_func++;
    if (indogr == 1)
    {
        if ((x >= -10) && (y >= 0) && (x + y <= 1))
        {
            return 0;
        }
        else
        {
            return 1;
        }
    }
    else
    {
        if (((x+3)*(x+3))/16+((y+4)*(y+4))/9<=1)
        {
            return 0;
        }
        else
        {
            return 1;
        }
    }
}
double NormVec(double x, double y)
{
    return sqrt(x * x + y * y);
}
std::vector<double> Simplex(double eps, double xn, double yn)
{
    std::vector<double> points;
    double  d, xc, b, yc, res, alph, bet, gam, xt[3], yt[3], xtres, ytres, f[3], fres, crit, delta;
    alph = 1;
    bet = 2;
    gam = 0.5;
    d = 0.5;
    delta = 1;
    xt[0] = xn;
    yt[0] = yn;
    xt[1] = xt[0] + l;
    yt[1] = yt[0] + l;
    xt[2] = xt[0] + l;
    yt[2] = yt[0];
    for (int i = 0; i < 3; i++)
    {
        f[i] = AimFunc(xt[i], yt[i]);
    }
    xc = (xt[1] + xt[2]) / 2;
    yc = (yt[1] + yt[2]) / 2;
    crit = 10;
    while (crit > eps * eps)
    {
        iter++;
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                if (f[j] < f[j + 1])
                {
                    res = f[j];
                    f[j] = f[j + 1];
                    f[j + 1] = res;
                    res = xt[j];
                    xt[j] = xt[j + 1];
                    xt[j + 1] = res;
                    res = yt[j];
                    yt[j] = yt[j + 1];
                    yt[j + 1] = res;
                }
            }
        }
        points.push_back(xt[0]);
        points.push_back(yt[0]);
        points.push_back(xt[1]);
        points.push_back(yt[1]);
        points.push_back(xt[2]);
        points.push_back(yt[2]);
        /*std::cout << xt[0] << ' ' << yt[0] << ' ' << f[0] << std::endl;
        std::cout << xt[1] << ' ' << yt[1] << ' ' << f[1] << std::endl;
        std::cout << xt[2] << ' ' << yt[2] << ' ' << f[2] << std::endl;
        std::cout<<std::endl;*/
        indu = 0;
        xtres = (xt[1] + xt[2]) / 2 + alph * ((xt[1] + xt[2]) / 2 - xt[0]);
        ytres = (yt[1] + yt[2]) / 2 + alph * ((yt[1] + yt[2]) / 2 - yt[0]);
        fres = AimFunc(xtres, ytres);
        if (CritOgr(xtres, ytres) == 0)
        {
            if (fres < f[0])
            {
                /*std::cout << iter << ' ' << 1 << std::endl;*/
                indu = 1;
                indc = 0;
                f[0] = fres;
                xt[0] = xtres;
                yt[0] = ytres;
            }
        }
        if (indu == 0)
        {
            xtres = (xt[2] + xt[0]) / 2 + alph * ((xt[2] + xt[0]) / 2 - xt[1]);
            ytres = (yt[2] + yt[0]) / 2 + alph * ((yt[2] + yt[0]) / 2 - yt[1]);
            if (CritOgr(xtres, ytres) == 0)
            {
                fres = AimFunc(xtres, ytres);
                if (fres < f[1])
                {
                    /*std::cout << iter << ' ' << 2 << std::endl;*/
                    indu = 1;
                    indc = 1;
                    f[1] = fres;
                    xt[1] = xtres;
                    yt[1] = ytres;
                }
            }
        }
        if (indu == 0)
        {
            xtres = (xt[0] + xt[1]) / 2 + alph * ((xt[0] + xt[1]) / 2 - xt[2]);
            ytres = (yt[0] + yt[1]) / 2 + alph * ((yt[0] + yt[1]) / 2 - yt[2]);
            fres = AimFunc(xtres, ytres);
            if (CritOgr(xtres, ytres) == 0)
            {
                if (fres < f[2])
                {
                    /*std::cout << iter << ' ' << 3 << std::endl;*/
                    indu = 1;
                    indc = 2;
                    f[2] = fres;
                    xt[2] = xtres;
                    yt[2] = ytres;
                }
            }
        }
        if (indu == 1)
        {

            xtres = (xt[(indc + 1) % 3] + xt[(indc + 2) % 3]) / 2 + bet * (xt[indc] - (xt[(indc + 1) % 3] + xt[(indc + 2) % 3]) / 2);
            ytres = (yt[(indc + 1) % 3] + yt[(indc + 2) % 3]) / 2 + bet * (yt[indc] - (yt[(indc + 1) % 3] + yt[(indc + 2) % 3]) / 2);
            fres = AimFunc(xtres, ytres);
            if (CritOgr(xtres, ytres) == 0)
            {
                if (f[indc] > fres)
                {
                    /*std::cout << iter << ' ' << "beta" << std::endl;*/
                    f[indc] = fres;
                    xt[indc] = xtres;
                    yt[indc] = ytres;
                }
            }
        }
        if (indu == 0)
        {
            xtres = (xt[1] + xt[2]) / 2 + gam * ((xt[1] + xt[2]) / 2 - xt[0]);
            ytres = (yt[1] + yt[2]) / 2 + gam * ((yt[1] + yt[2]) / 2 - yt[0]);
            if (CritOgr(xtres, ytres) == 0)
            {
                fres = AimFunc(xtres, ytres);
                if (f[0] > fres)
                {
                    /*std::cout << iter << ' ' << "gama1" <<fres << ' ' << f[0] << std::endl;*/
                    indu = 1;
                    f[0] = fres;
                    xt[0] = xtres;
                    yt[0] = ytres;
                }
            }
        }
        if (indu == 0)
        {
            xtres = (xt[1] + xt[2]) / 2 - gam * ((xt[1] + xt[2]) / 2 - xt[0]);
            ytres = (yt[1] + yt[2]) / 2 - gam * ((yt[1] + yt[2]) / 2 - yt[0]);
            fres = AimFunc(xtres, ytres);
            if (CritOgr(xtres, ytres) == 0)
            {
                if (f[0] > fres)
                {
                    /*std::cout << iter << ' ' << "gama2" <<fres<<' ' << f[0] << std::endl;*/
                    indu = 1;
                    f[0] = fres;
                    xt[0] = xtres;
                    yt[0] = ytres;
                }
            }
        }
        if (indu == 0)
        {
            /*xt[0] = (xt[1] + xt[2]) / 2 + gam * ((xt[1] + xt[2]) / 2 - xt[0]);
            yt[0] = (yt[1] + yt[2]) / 2 + gam * ((yt[1] + yt[2]) / 2 - yt[0]);*/
            delta = delta * d;
            xt[1] = xt[0] + delta * (xt[1] - xt[0]);
            yt[1] = yt[0] + delta * (yt[1] - yt[0]);
            xt[2] = xt[0] + delta * (xt[2] - xt[0]);
            yt[2] = yt[0] + delta * (yt[2] - yt[0]);
            for (int i = 1; i < 3; i++)
            {
                f[i] = AimFunc(xt[i], yt[i]);
            }
        }
        crit = 0.5 * fabs((xt[1] - xt[0]) * (yt[2] - yt[0]) - (xt[2] - xt[0]) * (yt[1] - yt[0]));
    }
    return points;
}
std::vector<double> Simplex2(double eps, double xn, double yn)
{
    std::vector<double> points;
    double  d, xc, b, yc, res, alph, bet, gam, xt[3], yt[3], xtres, ytres, f[3], fres, crit, delta;
    alph = 1;
    bet = 2;
    gam = 0.5;
    d = 0.5;
    delta = 1;
    xt[0] = xn;
    yt[0] = yn;
    xt[1] = xt[0] + l;
    yt[1] = yt[0] + l;
    xt[2] = xt[0] + l;
    yt[2] = yt[0];
    for (int i = 0; i < 3; i++)
    {
        f[i] = AimFunc(xt[i], yt[i]);
    }
    xc = (xt[1] + xt[2]) / 2;
    yc = (yt[1] + yt[2]) / 2;
    crit = 10;
    while (crit > eps * eps)
    {
        iter++;
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                if (f[j] < f[j + 1])
                {
                    res = f[j];
                    f[j] = f[j + 1];
                    f[j + 1] = res;
                    res = xt[j];
                    xt[j] = xt[j + 1];
                    xt[j + 1] = res;
                    res = yt[j];
                    yt[j] = yt[j + 1];
                    yt[j + 1] = res;
                }
            }
        }
        points.push_back(xt[0]);
        points.push_back(yt[0]);
        points.push_back(xt[1]);
        points.push_back(yt[1]);
        points.push_back(xt[2]);
        points.push_back(yt[2]);
        /*std::cout << xt[0] << ' ' << yt[0] << ' ' << f[0] << std::endl;
        std::cout << xt[1] << ' ' << yt[1] << ' ' << f[1] << std::endl;
        std::cout << xt[2] << ' ' << yt[2] << ' ' << f[2] << std::endl;
        std::cout<<std::endl;*/
        indu = 0;
        xtres = (xt[1] + xt[2]) / 2 + alph * ((xt[1] + xt[2]) / 2 - xt[0]);
        ytres = (yt[1] + yt[2]) / 2 + alph * ((yt[1] + yt[2]) / 2 - yt[0]);
        fres = AimFunc(xtres, ytres);
        if (fres < f[0])
        {
            /*std::cout << iter << ' ' << 1 << std::endl;*/
            indu = 1;
            indc = 0;
            f[0] = fres;
            xt[0] = xtres;
            yt[0] = ytres;
        }
        if (indu == 0)
        {
            xtres = (xt[2] + xt[0]) / 2 + alph * ((xt[2] + xt[0]) / 2 - xt[1]);
            ytres = (yt[2] + yt[0]) / 2 + alph * ((yt[2] + yt[0]) / 2 - yt[1]);
            fres = AimFunc(xtres, ytres);
            if (fres < f[1])
            {
                /*std::cout << iter << ' ' << 2 << std::endl;*/
                indu = 1;
                indc = 1;
                f[1] = fres;
                xt[1] = xtres;
                yt[1] = ytres;
            }
        }
        if (indu == 0)
        {
            xtres = (xt[0] + xt[1]) / 2 + alph * ((xt[0] + xt[1]) / 2 - xt[2]);
            ytres = (yt[0] + yt[1]) / 2 + alph * ((yt[0] + yt[1]) / 2 - yt[2]);
            fres = AimFunc(xtres, ytres);
            if (fres < f[2])
            {
                /*std::cout << iter << ' ' << 3 << std::endl;*/
                indu = 1;
                indc = 2;
                f[2] = fres;
                xt[2] = xtres;
                yt[2] = ytres;
            }
        }
        if (indu == 1)
        {

            xtres = (xt[(indc + 1) % 3] + xt[(indc + 2) % 3]) / 2 + bet * (xt[indc] - (xt[(indc + 1) % 3] + xt[(indc + 2) % 3]) / 2);
            ytres = (yt[(indc + 1) % 3] + yt[(indc + 2) % 3]) / 2 + bet * (yt[indc] - (yt[(indc + 1) % 3] + yt[(indc + 2) % 3]) / 2);
            fres = AimFunc(xtres, ytres);
            if (f[indc] > fres)
            {
                /*std::cout << iter << ' ' << "beta" << std::endl;*/
                f[indc] = fres;
                xt[indc] = xtres;
                yt[indc] = ytres;
            }
        }
        if (indu == 0)
        {
            xtres = (xt[1] + xt[2]) / 2 + gam * ((xt[1] + xt[2]) / 2 - xt[0]);
            ytres = (yt[1] + yt[2]) / 2 + gam * ((yt[1] + yt[2]) / 2 - yt[0]);
            fres = AimFunc(xtres, ytres);
            if (f[0] > fres)
            {
                /*std::cout << iter << ' ' << "gama1" <<fres << ' ' << f[0] << std::endl;*/
                indu = 1;
                f[0] = fres;
                xt[0] = xtres;
                yt[0] = ytres;
            }
        }
        if (indu == 0)
        {
            xtres = (xt[1] + xt[2]) / 2 - gam * ((xt[1] + xt[2]) / 2 - xt[0]);
            ytres = (yt[1] + yt[2]) / 2 - gam * ((yt[1] + yt[2]) / 2 - yt[0]);
            fres = AimFunc(xtres, ytres);
            if (f[0] > fres)
            {
                /*std::cout << iter << ' ' << "gama2" <<fres<<' ' << f[0] << std::endl;*/
                indu = 1;
                f[0] = fres;
                xt[0] = xtres;
                yt[0] = ytres;
            }
        }
        if (indu == 0)
        {
            /*xt[0] = (xt[1] + xt[2]) / 2 + gam * ((xt[1] + xt[2]) / 2 - xt[0]);
            yt[0] = (yt[1] + yt[2]) / 2 + gam * ((yt[1] + yt[2]) / 2 - yt[0]);*/
            delta = delta * d;
            xt[1] = xt[0] + delta * (xt[1] - xt[0]);
            yt[1] = yt[0] + delta * (yt[1] - yt[0]);
            xt[2] = xt[0] + delta * (xt[2] - xt[0]);
            yt[2] = yt[0] + delta * (yt[2] - yt[0]);
            for (int i = 1; i < 3; i++)
            {
                f[i] = AimFunc(xt[i], yt[i]);
            }
        }
        crit = 0.5 * fabs((xt[1] - xt[0]) * (yt[2] - yt[0]) - (xt[2] - xt[0]) * (yt[1] - yt[0]));
    }
    return points;
}
std::vector<double> RegSimplex(double eps, double xn, double yn)
{
    std::vector<double> points;
    double d, res, xt[3], crit, yt[3], xtres, ytres, f[3], fres, delta;
    d = 0.5;
    delta = 1;
    l = 1;
    indu = 0;
    xt[0] = xn;
    yt[0] = yn;
    xt[1] = xt[0] + l * ((sqrt(3) - 1) / (2 * sqrt(2)));
    yt[1] = yt[0] + l * ((sqrt(3) + 2 - 1) / (2 * sqrt(2)));
    xt[2] = xt[0] + l * ((sqrt(3) + 2 - 1) / (2 * sqrt(2)));
    yt[2] = yt[0] + l * ((sqrt(3) - 1) / (2 * sqrt(2)));
    for (int i = 0; i < 3; i++)
    {
        f[i] = AimFunc(xt[i], yt[i]);
    }
    crit = 10.;
    while (crit > eps * eps)
    {
        iter1++;
        for (int j = 1; j < 3; j++)
        {
            for (int i = 0; i < 2; i++)
            {
                if (f[i] < f[i + 1])
                {
                    res = f[i];
                    f[i] = f[i + 1];
                    f[i + 1] = res;
                    res = xt[i];
                    xt[i] = xt[i + 1];
                    xt[i + 1] = res;
                    res = yt[i];
                    yt[i] = yt[i + 1];
                    yt[i + 1] = res;
                }
            }
        }
        points.push_back(xt[0]);
        points.push_back(yt[0]);
        points.push_back(xt[1]);
        points.push_back(yt[1]);
        points.push_back(xt[2]);
        points.push_back(yt[2]);
        indu = 0;
        delta = 1;
        while (indu == 0)
        {
            xtres = 2 * ((xt[1] + xt[2]) / 2) - xt[0];
            ytres = 2 * ((yt[1] + yt[2]) / 2) - yt[0];
            fres = AimFunc(xtres, ytres);
            if (CritOgr(xtres, ytres) ==0)
            {
                if (fres < f[0])
                {
                    indu = 1;
                    f[0] = fres;
                    xt[0] = xtres;
                    yt[0] = ytres;
                }
            }
                 
            if (indu == 0)
            {
                xtres = 2 * ((xt[2] + xt[0]) / 2) - xt[1];
                ytres = 2 * ((yt[2] + yt[0]) / 2) - yt[1];
                fres = AimFunc(xtres, ytres);
                if (CritOgr(xtres, ytres) == 0)
                {
                    if (fres < f[1])
                    {
                        indu = 1;
                        f[1] = fres;
                        xt[1] = xtres;
                        yt[1] = ytres;
                    }
                }
            }
            if (indu == 0)
            {
                xtres = 2 * ((xt[0] + xt[1]) / 2) - xt[2];
                ytres = 2 * ((yt[0] + yt[1]) / 2) - yt[2];
                fres = AimFunc(xtres, ytres);
                if (CritOgr(xtres, ytres) == 0)
                {
                    if (fres < f[2])
                    {
                        indu = 1;
                        f[1] = fres;
                        xt[1] = xtres;
                        yt[1] = ytres;
                    }
                }
            }
            if (indu == 0)
            {
                delta = delta / 2;
                xt[1] = xt[0] + delta * (xt[1] - xt[0]);
                yt[1] = yt[0] + delta * (yt[1] - yt[0]);
                xt[2] = xt[0] + delta * (xt[2] - xt[0]);
                yt[2] = yt[0] + delta * (yt[2] - yt[0]);
                for (int i = 0; i < 3; i++)
                {
                    f[i] = AimFunc(xt[i], yt[i]);
                }
                indu = 1;
            }
        }
        crit = 0.5 * fabs((xt[1] - xt[0]) * (yt[2] - yt[0]) - (xt[2] - xt[0]) * (yt[1] - yt[0]));
    }
    /*std::cout << xt[0] << ' ' << yt[0] <<' ' <<f[0]<< std::endl;
   std::cout << xt[1] << ' ' << yt[1] << ' ' << f[1] << std::endl  ;
   std::cout << xt[2] << ' ' << yt[2] << ' ' << f[2] << std::endl ;*/
    return points;
}
std::vector<double> RegSimplex2(double eps, double xn, double yn)
{
    std::vector<double> points;
    double d, res, xt[3], crit, yt[3], xtres, ytres, f[3], fres, delta;
    d = 0.5;
    delta = 1;
    l = 1;
    indu = 0;
    xt[0] = xn;
    yt[0] = yn;
    xt[1] = xt[0] + l * ((sqrt(3) - 1) / (2 * sqrt(2)));
    yt[1] = yt[0] + l * ((sqrt(3) + 2 - 1) / (2 * sqrt(2)));
    xt[2] = xt[0] + l * ((sqrt(3) + 2 - 1) / (2 * sqrt(2)));
    yt[2] = yt[0] + l * ((sqrt(3) - 1) / (2 * sqrt(2)));
    for (int i = 0; i < 3; i++)
    {
        f[i] = AimFunc2(xt[i], yt[i]);
    }
    crit = 10.;
    while (crit > eps * eps)
    {
        iter1++;
        for (int j = 1; j < 3; j++)
        {
            for (int i = 0; i < 2; i++)
            {
                if (f[i] < f[i + 1])
                {
                    res = f[i];
                    f[i] = f[i + 1];
                    f[i + 1] = res;
                    res = xt[i];
                    xt[i] = xt[i + 1];
                    xt[i + 1] = res;
                    res = yt[i];
                    yt[i] = yt[i + 1];
                    yt[i + 1] = res;
                }
            }
        }
        points.push_back(xt[0]);
        points.push_back(yt[0]);
        points.push_back(xt[1]);
        points.push_back(yt[1]);
        points.push_back(xt[2]);
        points.push_back(yt[2]);
        indu = 0;
        delta = 1;
        while (indu == 0)
        {
            xtres = 2 * ((xt[1] + xt[2]) / 2) - xt[0];
            ytres = 2 * ((yt[1] + yt[2]) / 2) - yt[0];
            fres = AimFunc2(xtres, ytres);
            
                if (fres < f[0])
                {
                    indu = 1;
                    f[0] = fres;
                    xt[0] = xtres;
                    yt[0] = ytres;
                
            }

            if (indu == 0)
            {
                xtres = 2 * ((xt[2] + xt[0]) / 2) - xt[1];
                ytres = 2 * ((yt[2] + yt[0]) / 2) - yt[1];
                fres = AimFunc2(xtres, ytres);
                
                    if (fres < f[1])
                    {
                        indu = 1;
                        f[1] = fres;
                        xt[1] = xtres;
                        yt[1] = ytres;
                    }
                
            }
            if (indu == 0)
            {
                xtres = 2 * ((xt[0] + xt[1]) / 2) - xt[2];
                ytres = 2 * ((yt[0] + yt[1]) / 2) - yt[2];
                fres = AimFunc2(xtres, ytres);
                
                    if (fres < f[2])
                    {
                        indu = 1;
                        f[1] = fres;
                        xt[1] = xtres;
                        yt[1] = ytres;
                    }
                
            }
            if (indu == 0)
            {
                delta = delta / 2;
                xt[1] = xt[0] + delta * (xt[1] - xt[0]);
                yt[1] = yt[0] + delta * (yt[1] - yt[0]);
                xt[2] = xt[0] + delta * (xt[2] - xt[0]);
                yt[2] = yt[0] + delta * (yt[2] - yt[0]);
                for (int i = 0; i < 3; i++)
                {
                    f[i] = AimFunc2(xt[i], yt[i]);
                }
                indu = 1;
            }
        }
        crit = 0.5 * fabs((xt[1] - xt[0]) * (yt[2] - yt[0]) - (xt[2] - xt[0]) * (yt[1] - yt[0]));
    }
    /*std::cout << xt[0] << ' ' << yt[0] <<' ' <<f[0]<< std::endl;
   std::cout << xt[1] << ' ' << yt[1] << ' ' << f[1] << std::endl  ;
   std::cout << xt[2] << ' ' << yt[2] << ' ' << f[2] << std::endl ;*/
    return points;
}
std::vector<double> Vnutr(double eps, double xn, double yn)
{
    std::vector<double> points,p;
    double xt, yt, xtl, ytl;
    xtl = xn+1;
    ytl = yn+1;
    xt = xn;
    yt = yn;
    points.push_back(xt);
    points.push_back(yt);
    p = RegSimplex(eps, xn, yn);
    xt = p[p.size() - 2];
    yt = p[p.size() - 1];
    points.push_back(xt);
    points.push_back(yt);
    while (fabs(AimFunc(xt,yt)- AimFunc(xtl, ytl))>eps/10)
    {
        iter++;
        r = r * betta;
        xtl = xt;
        ytl = yt;
        p = RegSimplex(eps/100, xn, yn);
        xt = (p[p.size() - 2]+ p[p.size() - 4]+ p[p.size() - 6])/3;
        yt = (p[p.size() - 1] + p[p.size() - 3] + p[p.size() - 5]) / 3;
        points.push_back(xt);
        points.push_back(yt);
        if (iter > 50)
        {
            break;
        }
    }
    return points;
}
std::vector<double> Vnesh(double eps, double xn, double yn)
{
    std::vector<double> points, p;
    double xt, yt, xtl, ytl;
    xtl = xn + 1;
    ytl = yn + 1;
    xt = xn;
    yt = yn;
    points.push_back(xt);
    points.push_back(yt);
    p = RegSimplex2(eps, xn, yn);
    xt = p[p.size() - 2];
    yt = p[p.size() - 1];
    points.push_back(xt);
    points.push_back(yt);
    while (fabs(AimFunc2(xt, yt) - AimFunc2(xtl, ytl)) > eps )
    {
        iter++;
        r = r +50;
        xtl = xt;
        ytl = yt;
        p = RegSimplex2(eps , xn, yn);
        xt = (p[p.size() - 2] + p[p.size() - 4] + p[p.size() - 6]) / 3;
        yt = (p[p.size() - 1] + p[p.size() - 3] + p[p.size() - 5]) / 3;
        points.push_back(xt);
        points.push_back(yt);
        if (iter > 1000)
        {
            break;
        }
        std::cout << AimFunc2(xt, yt) << ' ' << AimFunc2(xtl, ytl) << std::endl;
    }
    return points;
}
int main()
{
    p = 3;
    l = 1;
    const double eps[2] = { 0.01,0.00001 };;
    std::vector<double> points1, points2;
    std::ofstream file;
    std::string fil, end;
    std::cout << std::setprecision(7);
    file << std::setprecision(7);
    indogr = 1;
    betta = 0.5;
    double x0[2][2] = { { -7,-9},{-7,0} };
    double y0[2][2] = { {4,2},{4,-3} };
    //Constobn = 4;  // Константа обновления
    indogr =1; //Ограничения
    indf = 3;   //функция
    int i = 2;   //точность
    int j =1; //точка
    fil = "Vnutreps";
    end = std::to_string(i);
    end += "f";
    end += std::to_string(indf);
    end += "Ogr";
    end += std::to_string(indogr);
    end += "p";
    end += std::to_string(j);
    end += ".txt";
    std::cout << std::setprecision(7);
    file << std::setprecision(7);
    file.open(fil + end);
    culc_func = 0;
    iter = 0;
    r = 32;
    l = 2;
    //r = 32;
    points1 = Vnutr(eps[i - 1], x0[indogr - 1][j-1], y0[indogr - 1][j-1]);
    file << iter << ' ' << culc_func << std::endl;
    file << std::setprecision(7);
    for (int i = 0; i < points1.size() - 1; i += 2)
    {
        file << points1[i] << ' ' << points1[i + 1] << std::endl;
    }
    file.close();
    std::cout<<"Vnutr: " << iter << ' ' << culc_func << std::endl;
    std::cout << "Point: " << points1[points1.size()-2] << ' ' << points1[points1.size() - 1] << std::endl;
    std::cout << std::endl;
    fil = "Vnesheps";
    end = std::to_string(i);
    end += "f";
    end += std::to_string(indf);
    end += "Ogr";
    end += std::to_string(indogr);
    end += "p";
    end += std::to_string(j);
    end += ".txt";
    std::cout << std::setprecision(7);
    file << std::setprecision(7);
    file.open(fil + end);
    culc_func = 0;
    iter = 0;
    r = 2;
    l = 2;
    //r = 32;
    points1 = Vnesh(eps[i - 1], x0[indogr - 1][j - 1], y0[indogr - 1][j - 1]);
    file << iter << ' ' << culc_func << std::endl;
    file << std::setprecision(7);
    for (int i = 0; i < points1.size() - 1; i += 2)
    {
        file << points1[i] << ' ' << points1[i + 1] << std::endl;
    }
    file.close();
    std::cout << "Vnesh: " << iter << ' ' << culc_func << std::endl;
    std::cout << "Point: " << points1[points1.size() - 2] << ' ' << points1[points1.size() - 1] << std::endl;
    std::cout << std::endl;

    
    std::cout << "Hello World!\n";
    return 0;
}

