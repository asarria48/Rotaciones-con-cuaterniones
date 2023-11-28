#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <stdio.h>
using namespace std;

double F1(double, double);
double F2(double, double);
double F3(double, double);
double ff0(double, double, double, double, double, double);
double ff1(double, double, double, double, double, double, double);
double ff2(double, double, double, double, double, double, double);
double ff3(double, double, double, double, double, double, double);
double QQ(double, double, double, double);

//Constantes

const double m = 2;
const double M = 4;
const double R = 2;

const double I1 = 2*m*R*R;
const double I2 = 2*M*R*R;
const double I3 = 2*(m+M)*R*R;

//const double I1 = 1;
//const double I2 = 2;
//const double I3 = 3;

const double pi = 3.14159265;

int main() {
  //Abrimos los archivos dat
  ofstream outfilew;
  ofstream outfileq;

  outfilew.open("Datos-w.dat");
  outfileq.open("Datos-q.dat");

  //Datos iniciales de las velocidades angulares y el cuaternión
  double w10 = 0.3;
  double w20 = 30;
  double w30 = 0.35;
  double q00 = 1;
  double q10 = 0;
  double q20 = 0;
  double q30 = 0;

  //Precisión
  double h=0.01;
  //Tiempos
  double ti=0.0;
  double tf=2.0;
  //Número de pasos
  int NPasos=(tf-ti)/h;
  //Arreglos de w y q
  double t[NPasos+1]={};
  double w1[NPasos+1]={};
  double w2[NPasos+1]={};
  double w3[NPasos+1]={};
  double q0[NPasos+1]={};
  double q1[NPasos+1]={};
  double q2[NPasos+1]={};
  double q3[NPasos+1]={};

  //Inicializamos y exportamos los primeros valores de los arreglos
  t[0]=ti;
  w1[0]=w10;
  w2[0]=w20;
  w3[0]=w30;
  q0[0]=q00;
  q1[0]=q10;
  q2[0]=q20;
  q3[0]=q30;
  outfilew << t[0] << "," << w1[0] << "," << w2[0]  << "," << w3[0] <<  endl;
  outfileq << t[0] << "," << q0[0] << "," << q1[0] << "," << q2[0] << "," << q3[0] <<  endl;

  //Calculamos los valores respectivos de cada arreglo con RK4
  for (int i=0; i<NPasos; i++){
    //tiempo
    t[i+1]=t[i]+h;

    //K1
    double k1w1 = h * F1(w2[i], w3[i]);
    double k1w2 = h * F2(w1[i], w3[i]);
    double k1w3 = h * F3(w1[i], w2[i]);
    double k1q0 = h * ff0(w1[i], w2[i], w3[i], q1[i], q2[i], q3[i]);
    double k1q1 = h * ff1(w1[i], w2[i], w3[i], q0[i], q1[i], q2[i], q3[i]);
    double k1q2 = h * ff2(w1[i], w2[i], w3[i], q0[i], q1[i], q2[i], q3[i]);
    double k1q3 = h * ff3(w1[i], w2[i], w3[i], q0[i], q1[i], q2[i], q3[i]);

    //K2 (tener en cuenta que a cada variable se suma dentro de la función respectiva el k anterior correspondiente a esa misma variable)
    double k2w1 = h * F1(w2[i] + 0.5 * k1w2, w3[i] + 0.5 * k1w3);
    double k2w2 = h * F2(w1[i] + 0.5 * k1w1, w3[i] + 0.5 * k1w3);
    double k2w3 = h * F3(w1[i] + 0.5 * k1w1, w2[i] + 0.5 * k1w2);
    double k2q0 = h * ff0(w1[i] + 0.5 * k1w1, w2[i] + 0.5 * k1w2, w3[i] + 0.5 * k1w3, q1[i] + 0.5 * k1q1, q2[i] + 0.5 * k1q2, q3[i] + 0.5 * k1q3);
    double k2q1 = h * ff1(w1[i] + 0.5 * k1w1, w2[i] + 0.5 * k1w2, w3[i] + 0.5 * k1w3, q0[i] + 0.5 * k1q0, q1[i] + 0.5 * k1q1, q2[i] + 0.5 * k1q2, q3[i] + 0.5 * k1q3);
    double k2q2 = h * ff2(w1[i] + 0.5 * k1w1, w2[i] + 0.5 * k1w2, w3[i] + 0.5 * k1w3, q0[i] + 0.5 * k1q0, q1[i] + 0.5 * k1q1, q2[i] + 0.5 * k1q2, q3[i] + 0.5 * k1q3);
    double k2q3 = h * ff3(w1[i] + 0.5 * k1w1, w2[i] + 0.5 * k1w2, w3[i] + 0.5 * k1w3, q0[i] + 0.5 * k1q0, q1[i] + 0.5 * k1q1, q2[i] + 0.5 * k1q2, q3[i] + 0.5 * k1q3);

    //K3
    double k3w1 = h * F1(w2[i] + 0.5 * k2w2, w3[i] + 0.5 * k2w3);
    double k3w2 = h * F2(w1[i] + 0.5 * k2w1, w3[i] + 0.5 * k2w3);
    double k3w3 = h * F3(w1[i] + 0.5 * k2w1, w2[i] + 0.5 * k2w2);
    double k3q0 = h * ff0(w1[i] + 0.5 * k2w1, w2[i] + 0.5 * k2w2, w3[i] + 0.5 * k2w3, q1[i] + 0.5 * k2q1, q2[i] + 0.5 * k2q2, q3[i] + 0.5 * k2q3);
    double k3q1 = h * ff1(w1[i] + 0.5 * k2w1, w2[i] + 0.5 * k2w2, w3[i] + 0.5 * k2w3, q0[i] + 0.5 * k2q0, q1[i] + 0.5 * k2q1, q2[i] + 0.5 * k2q2, q3[i] + 0.5 * k2q3);
    double k3q2 = h * ff2(w1[i] + 0.5 * k2w1, w2[i] + 0.5 * k2w2, w3[i] + 0.5 * k2w3, q0[i] + 0.5 * k2q0, q1[i] + 0.5 * k2q1, q2[i] + 0.5 * k2q2, q3[i] + 0.5 * k2q3);
    double k3q3 = h * ff3(w1[i] + 0.5 * k2w1, w2[i] + 0.5 * k2w2, w3[i] + 0.5 * k2w3, q0[i] + 0.5 * k2q0, q1[i] + 0.5 * k2q1, q2[i] + 0.5 * k2q2, q3[i] + 0.5 * k2q3);

    //K4
    double k4w1 = h * F1(w2[i] + k3w2, w3[i] + k3w3);
    double k4w2 = h * F2(w1[i] + k3w1, w3[i] + k3w3);
    double k4w3 = h * F3(w1[i] + k3w1, w2[i] + k3w2);
    double k4q0 = h * ff0(w1[i] + k3w1, w2[i] + k3w2, w3[i] + k3w3, q1[i] + k3q1, q2[i] + k3q2, q3[i] + k3q3);
    double k4q1 = h * ff1(w1[i] + k3w1, w2[i] + k3w2, w3[i] + k3w3, q0[i] + k3q0, q1[i] + k3q1, q2[i] + k3q2, q3[i] + k3q3);
    double k4q2 = h * ff2(w1[i] + k3w1, w2[i] + k3w2, w3[i] + k3w3, q0[i] + k3q0, q1[i] + k3q1, q2[i] + k3q2, q3[i] + k3q3);
    double k4q3 = h * ff3(w1[i] + k3w1, w2[i] + k3w2, w3[i] + k3w3, q0[i] + k3q0, q1[i] + k3q1, q2[i] + k3q2, q3[i] + k3q3);

    // Calculamos y exportamos w1, w2, w3
    w1[i+1]=w1[i]+(k1w1+2*k2w1+2*k3w1+k4w1)/6;
    w2[i+1]=w2[i]+(k1w2+2*k2w2+2*k3w2+k4w2)/6;
    w3[i+1]=w3[i]+(k1w3+2*k2w3+2*k3w3+k4w3)/6;
    outfilew << t[i+1] << "," << w1[i+1] << "," << w2[i+1] << "," << w3[i+1] << endl;

    //calculamos q0, q1, q2, q3
    q0[i+1]=q0[i]+(k1q0+2*k2q0+2*k3q0+k4q0)/6;
    q1[i+1]=q1[i]+(k1q1+2*k2q1+2*k3q1+k4q1)/6;
    q2[i+1]=q2[i]+(k1q2+2*k2q2+2*k3q2+k4q2)/6;
    q3[i+1]=q3[i]+(k1q3+2*k2q3+2*k3q3+k4q3)/6;

    //normalizamos el cuaternión
    double Q = QQ(q0[i+1], q1[i+1], q2[i+1], q3[i+1]);
    q0[i+1]=q0[i+1]/Q;
    q1[i+1]=q1[i+1]/Q;
    q2[i+1]=q2[i+1]/Q;
    q3[i+1]=q3[i+1]/Q;

    //Exportamos q0, q1, q2, q3
    outfileq << t[i+1] << "," << q0[i+1] << "," << q1[i+1] << "," << q2[i+1] << "," << q3[i+1] << endl;

  }

  //Cerramos los archivos dat
  outfilew.close();
  outfileq.close();

  return 0; 
}

//Funciones para calcular w
double F1(double w2, double w3){

  return -w2*w3*(I3 - I2)/I1;
}

double F2(double w1, double w3){

  return -w1*w3*(I1 - I3)/I2;
}

double F3(double w1, double w2){

  return -w1*w2*(I2 - I1)/I3;
}

//Funciones para calcular q
double ff0(double w1, double w2, double w3, double q1, double q2, double q3){

  return -0.5*(w1*q1 + w2*q2 + w3*q3);
}

double ff1(double w1, double w2, double w3, double q0, double q1, double q2, double q3){

  return 0.5*(q0*w1 + w2*q3 - w3*q2);
}

double ff2(double w1, double w2, double w3, double q0, double q1, double q2, double q3){

  return 0.5*(q0*w2 - w1*q3 + w3*q1);
}

double ff3(double w1, double w2, double w3, double q0, double q1, double q2, double q3){

  return 0.5*(q0*w3 + w1*q2 - w2*q1);
}

//Función de normalización de q
double QQ(double q0, double q1, double q2, double q3){

  return std::sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);

}