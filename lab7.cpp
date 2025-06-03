#include<iostream>
#include <fstream>
#include<cmath>
#include<math.h>
#include<vector>
#include <iomanip>
#include<string>
int iter = 0;
int culc_func = 0;
const double alpha1 = 1.0;
const double alpha2 = 5.0;
int Constobn = 2;
int indf = 3;
int indu = 0;
int indc = 0;
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
    if (indf == 4)
    {
        return (x * x + y - 11) * (x * x + y - 11) + (x + y * y - 7) * (x + y * y - 7);
    }
}
double NormVec(double x, double y)
{
    return sqrt(x * x + y * y);
}
double CritSimplex(double xc, double yc,double* xt, double* yt)
{
    double res,b;
    res = 0;
    for (int i = 0; i < 3; i++)
    {
        b = AimFunc(xt[i], yt[i]) - AimFunc(xc, yc);
        res = res + b * b;
    }
    return (1/3)*res;
}
std::vector<double> RegSimplex(double eps, double xn, double yn)
{
    std::vector<double> points;
    double d, res,xt[3],crit, yt[3], xtres, ytres, f[3],fres, delta;
    d = 0.5;
    delta = 1;
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
    while (crit > eps*eps)
    {
        iter++;
        for (int j = 1; j < 3; j++)
        {
            for (int i = 0; i < 2; i++)
            {
                if (f[i] < f[i+1])
                {
                    res = f[i];
                    f[i] = f[i+1];
                    f[i+1] = res;
                    res = xt[i];
                    xt[i] = xt[i+1];
                    xt[i+1] = res;
                    res = yt[i];
                    yt[i] = yt[i+1];
                    yt[i+1] = res;
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
                fres = AimFunc(xtres, ytres);
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
                fres = AimFunc(xtres, ytres);
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
                delta = delta /2;
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
std::vector<double> Simplex(double eps, double xn, double yn)
{
    std::vector<double> points;
    double  d,xc,b,yc ,res,alph,bet,gam, xt[3], yt[3], xtres, ytres, f[3], fres,crit, delta;
    alph = 1;
    bet = 2;
    gam = 0.5;
    d = 0.5;
    delta = 1;
    xt[0] = xn;
    yt[0] = yn;
    xt[1] = xt[0] +l;
    yt[1] = yt[0] + l ;
    xt[2] = xt[0] + l;
    yt[2] = yt[0] ;
    for (int i = 0; i < 3; i++)
    {
        f[i] = AimFunc(xt[i], yt[i]);
    }
    xc = (xt[1] + xt[2]) / 2;
    yc = (yt[1] + yt[2]) / 2;
    crit = 10;
    while (crit> eps * eps)
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
            
            xtres = (xt[(indc+1)%3] + xt[(indc + 2) % 3]) / 2 + bet * (xt[indc] - (xt[(indc + 1) % 3] + xt[(indc + 2) % 3]) / 2);
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
        if(indu==0)
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
int main()
{
    const double eps[2] = { 0.01,0.00001 };
    std::vector<double> points1, points2;
    std::ofstream file;
    std::string fil, end;
    double x0[2] = { 2,10 };
    double y0[2] = { -5,-10 };
    p = 3;
    l = 1;
    //Constobn = 4;  // Константа обновления
    indf = 3;   //функция
    int i =1;   //точность
    int j = 2; //точка
    fil = "RegSimplexeps";
    end = std::to_string(i);
    end += "f";
    end += std::to_string(indf);
    end += "p";
    end += std::to_string(j);
    end += ".txt";
    std::cout << std::setprecision(15);
    file << std::setprecision(15);
    file.open(fil + end);
    iter = 0;
    culc_func = 0;
    points1 = RegSimplex(eps[i - 1], x0[j - 1], y0[j - 1]);
    file << iter << ' ' << culc_func << std::endl;
    file << std::setprecision(7);
    for (int i = 0; i < points1.size() - 1; i += 6)
    {
        file << points1[i] << ' ' << points1[i + 1]<< ' ' << points1[i + 2]<< ' ' << points1[i + 3]<< ' ' << points1[i + 4]<< ' ' << points1[i + 5] << std::endl;
    }
    file.close();
    std::cout << "RegSimplex:" << std::endl;
    std::cout << iter << ' ' << culc_func << std::endl;
    std::cout << std::endl;
    fil = "Simplexeps";
    end = std::to_string(i);
    end += "f";
    end += std::to_string(indf);
    end += "p";
    end += std::to_string(j);
    end += ".txt";
    std::cout << std::setprecision(15);
    file << std::setprecision(15);
    file.open(fil + end);
    iter = 0;
    culc_func = 0;
    points1 = Simplex(eps[i - 1], x0[j - 1], y0[j - 1]);
    file << iter << ' ' << culc_func << std::endl;
    file << std::setprecision(7);
    for (int i = 0; i < points1.size() - 1; i += 6)
    {
        file << points1[i] << ' ' << points1[i + 1] << ' ' << points1[i + 2] << ' ' << points1[i + 3] << ' ' << points1[i + 4] << ' ' << points1[i + 5] << std::endl;
    }
    file.close();
    std::cout << "Simplex:" << std::endl;
    std::cout << iter << ' ' << culc_func << std::endl;
    std::cout << std::endl;
    std::cout << "Hello World!\n";
    return 0;
}
