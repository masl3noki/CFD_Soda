#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "C:\MinGW\include\Eigen\Dense"


using namespace std;
using namespace Eigen;



Vector3d
from_U_to_F(Vector3d U ,double delta_x, double delta_t, double gamma){
    double rho = U[0];//считаем потоки через вектор U
    double u = U[1]/U[0];
    double E = U[2];
    double Gamma = gamma;
    Vector3d F;
    F[0] = U[1];
    F[1] = (E-0.5*rho*u*u)*(Gamma - 1.)+ rho*u*u;
    F[2] = u*(E + (E-0.5*rho*u*u)*(Gamma-1.));
    return F;
}
double
signum(double number){
    if(number>0.)
        return 1.;
    if(number<0.)
        return -1.;
    if(number == 0.)
        return 0;

    return 0;

}
Vector3d
found_alpha(Vector3d U_m_plus_one, Vector3d U_m, Vector3d U_m_minus_one,double delta_x){ //минмод
    Vector3d alpha;
    Vector3d u_m_plus_one = U_m_plus_one;// формулы 2.7.30 и 2.7.31
    Vector3d u_m_minus_one = U_m_minus_one;
    Vector3d u_m = U_m;
    double Delta_x = delta_x;
    for(int i = 0;i<3;i++){
        alpha[i]=0.5*((signum((u_m_plus_one[i]-u_m[i])/Delta_x))+signum((u_m[i]-u_m_minus_one[i])/Delta_x))*min(abs((u_m_plus_one[i]-u_m[i])/Delta_x),abs((u_m[i]-u_m_minus_one[i])/Delta_x));
    }
    //alpha[0]= max(0.5*(u_m_plus_one[0] - u_m[0] - abs(u_m_plus_one[0] - u_m[0])),0.5*(u_m[0] - u_m_minus_one[0] - abs(u_m[0] - u_m_minus_one[0])))+min(0.5*(u_m_plus_one[0] - u_m[0] + abs(u_m_plus_one[0] - u_m[0])),0.5*(u_m[0] - u_m_minus_one[0] + abs(u_m[0] - u_m_minus_one[0])));
    //alpha[1]= max(0.5*(u_m_plus_one[1] - u_m[1] - abs(u_m_plus_one[1] - u_m[1])),0.5*(u_m[1] - u_m_minus_one[1] - abs(u_m[1] - u_m_minus_one[1])))+min(0.5*(u_m_plus_one[1] - u_m[1] + abs(u_m_plus_one[1] - u_m[1])),0.5*(u_m[1] - u_m_minus_one[1] + abs(u_m[1] - u_m_minus_one[1])));
    //alpha[2]= max(0.5*(u_m_plus_one[2] - u_m[2] - abs(u_m_plus_one[2] - u_m[2])),0.5*(u_m[2] - u_m_minus_one[2] - abs(u_m[2] - u_m_minus_one[2])))+min(0.5*(u_m_plus_one[2] - u_m[2] + abs(u_m_plus_one[2] - u_m[2])),0.5*(u_m[2] - u_m_minus_one[2] + abs(u_m[2] - u_m_minus_one[2])));

    return alpha;
}

Vector3d
TVD_Godunov(Vector3d U_m_plus_one, Vector3d U_m, Vector3d U_m_minus_one, Vector3d U_m_minus_two,Vector3d U_m_plus_two, double delta_x, double delta_t, double gamma, Matrix3d A, Matrix3d Aabs ){
    Vector3d U_with_hat;
    Matrix3d a = A; //здесь просто переписываем все полученный переменные и исользуем формулы 2.5.5-2.5.10, а также для итоговой формулы использовал формулы из лекции на стр.50
    Matrix3d aabs = Aabs;
    Vector3d U_center;
    Vector3d u_m_plus_one = U_m_plus_one;
    //cout<<"UUUU==="<<u_m_plus_one<<endl;
    Vector3d u_m_minus_one = U_m_minus_one;
    Vector3d u_m_minus_two = U_m_minus_two;
    Vector3d u_m_plus_two = U_m_plus_two;
    Vector3d u_m = U_m;
    double Delta_x = delta_x;
    double Delta_t = delta_t;
    double Gamma = gamma;

    Vector3d U_l1 = u_m_minus_one+0.5*Delta_x*found_alpha(u_m,u_m_minus_one,u_m_minus_two,Delta_x);// считаем значения функций на границах
    Vector3d U_r1 = u_m-0.5*Delta_x*found_alpha(u_m_plus_one,u_m,u_m_minus_one,Delta_x);
    Vector3d U_l2 = u_m+0.5*Delta_x*found_alpha(u_m_plus_one,u_m,u_m_minus_one,Delta_x);
    Vector3d U_r2 = u_m_plus_one-0.5*Delta_x*found_alpha(u_m_plus_two,u_m_plus_one,u_m,Delta_x);

    U_with_hat = u_m+Delta_t/Delta_x*(0.5*(a*U_l1+a*U_r1+aabs*(U_l1-U_r1))-0.5*(a*U_l2+a*U_r2+aabs*(U_l2-U_r2))); //промежуточное значение функции 2.5.7
    //U_with_hat = u_m+ Delta_t/Delta_x*(from_U_to_F((u_m - 0.5*Delta_x*found_alpha(u_m_plus_one, u_m, u_m_minus_one,Delta_x)),Delta_x,Delta_t,Gamma) - from_U_to_F((u_m + 0.5*Delta_x*found_alpha(u_m_plus_one, u_m,u_m_minus_one,Delta_x)),Delta_x,Delta_t,Gamma));
    U_center = 0.5*(U_with_hat + u_m); //2.5.8
    //cout<<"u="<<0.5*Delta_x*found_alpha(u_m_plus_one, u_m, u_m_minus_one,Delta_x)<<endl;
    //cout<<"u_m+...="<<from_U_to_F((U_center - 0.5*Delta_x*found_alpha( u_m_plus_one, u_m, u_m_minus_one)),Delta_x,Delta_t,Gamma)<<endl;

    U_l1 = u_m_minus_one+0.5*Delta_x*found_alpha(u_m,u_m_minus_one,u_m_minus_two,Delta_x); 
    U_r1 = U_center-0.5*Delta_x*found_alpha(u_m_plus_one,u_m,u_m_minus_one,Delta_x);
    U_l2 = U_center+0.5*Delta_x*found_alpha(u_m_plus_one,u_m,u_m_minus_one,Delta_x);
    U_r2 = u_m_plus_one-0.5*Delta_x*found_alpha(u_m_plus_two,u_m_plus_one,u_m,Delta_x);

    return u_m+Delta_t/Delta_x*(0.5*(from_U_to_F(U_l1,Delta_x, Delta_t, Gamma)+from_U_to_F (U_r1,Delta_x, Delta_t, Gamma)+aabs*(U_l1-U_r1))-0.5*(from_U_to_F(U_l2,Delta_x, Delta_t, Gamma)+from_U_to_F(U_r2,Delta_x, Delta_t, Gamma)+aabs*(U_l2-U_r2)));
    //return u_m+Delta_t/Delta_x*(0.5*(from_U_to_F(u_m_minus_one,Delta_x, Delta_t, Gamma)+from_U_to_F(u_m,Delta_x, Delta_t, Gamma)+aabs*(u_m_minus_one-u_m))-0.5*(from_U_to_F(u_m,Delta_x, Delta_t, Gamma)+from_U_to_F(u_m_plus_one,Delta_x, Delta_t, Gamma)+aabs*(u_m-u_m_plus_one)));

}

// не используется
Vector3d
KIR(Vector3d Wl_1n, Vector3d Wln, Vector3d Wl__1n, Matrix3d OmegaT, Matrix3d A, Matrix3d H, double h, double tau) {
    return Wln - tau * A * (Wl__1n - Wl_1n) / (2 * h) +
           tau * (OmegaT.inverse() * H * OmegaT) * (Wl__1n - 2 * Wln + Wl_1n) / (2 * h);
}


int main() {

    //Неизвестные
    double u;
    double e;
    //

    double L = 1.; //длина области
    double h = L / 100.;//шаг
    const int Nx = 200;//количество делений нашей области, если менять, то надо поменять и выше
    double tau = 0.000001;// шаг по времени
    double g = 1.4;// это гамма
    double Max_uc = 0.;//максимальная скорость в области, нужна для определения длины шага дальше
    double Co = 0.005;
    vector<double> lol;
    double c;// speed of sound
    double t = 0.;// текущее время рассчета
    double T = 0.15;// это время всего рассчета
    int j = 0;
    Eigen::array<Vector3d, Nx> w0;//вектор решения на данном шаге
    Eigen::array<Vector3d, Nx> w1;// вектор решения на следующем шаге

    // нулевой слой, начальные условия, плотность, скорость, давление, слева и справа соответственно; НУ
    for (int i = 0; i < Nx / 2; i++) {
        w0[i][0] = 1.0;
        w0[i][1] = -2.0;
        w0[i][2] = 0.4 / (g - 1.);
    }
    
    for (int i = Nx / 2; i < Nx; i++) {
        w0[i][0] = 1.0;
        w0[i][1] = 2.0;
        w0[i][2] = 0.4 / (g - 1.);
    }
    Matrix3d A;
    Matrix3d Aabs;
    Matrix3d H;
    Matrix3d Habs;
    Matrix3d OmegaT;

    while (t <= T) {
        for (int i = 0; i < Nx; i++) {
            u = w0[i][1] / w0[i][0];
            e = w0[i][2] / w0[i][0]-0.5*u*u;
            c = sqrt((g - 1.) * g * e);//ищем скорость для сравнения с максимальной
            //if (max(max(abs(u + c), u), abs(u - c)) > Max_uc) {
               // Max_uc = max(max(abs(u + c), u), abs(u - c));
           // }
            OmegaT << -u * c, c, g - 1.,
                    -c * c, 0., g - 1.,
                    u * c, -c, g - 1.;
            A << 0., 1., 0.,
                    -u * u, 2. * u, g - 1.,
                    -g * u * e, g * e, u;
            //Aabs << 0, 1, 0,
             //       abs(-u * u),abs( 2 * u), g - 1,
              //      abs(-g * u * e), abs(g * e),abs( u);
            
            Habs << abs(u + c), 0., 0.,
                    0., abs(u), 0.,
                    0., 0., abs(u - c);
            H << u + c, 0., 0.,
                    0., u, 0.,
                    0., 0., u - c;

            Aabs = OmegaT.inverse()*Habs*OmegaT;
            if ((i > 1) && (i < Nx - 1)) {
                w1[i] = TVD_Godunov(w0[i + 1], w0[i], w0[i - 1],w0[i-2],w0[i+2], h, tau,g,OmegaT.inverse()*H*OmegaT,Aabs);//рассчет в области
            } else if (i == 0) {
                w1[i] = TVD_Godunov(w0[i+1], w0[i], w0[i],w0[i],w0[i+2], h, tau, g,OmegaT.inverse()*H*OmegaT,Aabs);// расчет на левой границе и ниже на правой
            } else if (i == Nx - 1) {
                w1[i] = TVD_Godunov(w0[i], w0[i], w0[i-1],w0[i-2],w0[i], h, tau,g,OmegaT.inverse()*H*OmegaT,Aabs);
            } else if (i==1){
                w1[i] = TVD_Godunov(w0[i+1], w0[i], w0[i-1],w0[i-1],w0[i+2], h, tau,g,OmegaT.inverse()*H*OmegaT,Aabs);
            }
            else if (i == Nx-2){
                w1[i] = TVD_Godunov(w0[i+1], w0[i], w0[i-1],w0[i-2],w0[i+1], h, tau,g,OmegaT.inverse()*H*OmegaT,Aabs);
            }
            
        }
        for (int i = 0; i < Nx; i++) {
            w0[i] = w1[i];//переопределили
        }
        lol.push_back(tau);
        cout <<"t="<< t<<"|||tau="<<tau << endl;
        t += tau;
        //tau = Co * h / Max_uc;//пересчитали шаг
    }
    Eigen::array<double, Nx> P;
    Eigen::array<double, Nx> U;
    Eigen::array<double, Nx> Ro;
    Eigen::array<double, Nx> E;
    for (int i = 0; i < Nx; i++) {
        P[i] = w1[i][2] * (g - 1);
        U[i] = w1[i][1] / w1[i][0];
        Ro[i] = w1[i][0];
        E[i] = w1[i][2] / w1[i][0];
    }
    ofstream out("P1.txt");
    for (int x = 0; x < Nx; x++) {
        out << P[x] / 100000 << ' ';
    }
    ofstream out1("U.txt");
    for (int x = 0; x < Nx; x++) {
        out1 << U[x] << ' ';
    }
    ofstream out2("Ro.txt");
    for (int x = 0; x < Nx; x++) {
        out2 << Ro[x] << ' ';
    }
    ofstream out3("E.txt");
    for (int x = 0; x < Nx; x++) {
        out3 << E[x] << ' ';
    }
    ofstream out4("T.txt");
    for (int x = 0; x < lol.size(); x++) {
        out4 << lol[x]*1000000 << ' ';
    }
}




















