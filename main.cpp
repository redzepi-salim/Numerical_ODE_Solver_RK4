#include <iostream>
#include <vector>
#include <string>
#include <functional>
#include "pugixml.hpp"
#include "matplotlibcpp.h"

using namespace std;
using namespace pugi;
namespace plt = matplotlibcpp;

plt::Plot plot("plot", "b-");

const double PI = 4 * atan(1);

//Funkcija koja pretvara string u decimalni broj
double Broj(std::string x) { return std::stod(x); }

// Funkcija koja trazi koeficijente dif. jed. iz stringa
std::vector<double> PronadjiKoeficijente(std::string s) {
    std::vector<double> koeficijenti;
    for (int i = 0; i < s.size(); i++) {
        int pocetak, kraj;
        bool predznak = false;
        std::string temp;
        if (s.at(i) >= '0' && s.at(i) <= '9') {
            pocetak = i;
            if (s.at(i - 1) == '-' || s.at(i - 2) == '-') predznak = true;
            i++;
            while ((s.at(i) >= '0' && s.at(i) <= '9') || s.at(i) == '.') i++;
            kraj = i;
            while (pocetak != (kraj + 1)) {
                temp.push_back(s.at(pocetak));
                pocetak++;
            }
            koeficijenti.push_back(Broj(temp));
            if (predznak) {
                koeficijenti.at(koeficijenti.size() - 1) = -koeficijenti.at(koeficijenti.size() - 1);
                predznak = false;
            }
        }
    }
    return koeficijenti;
}

//Funkcija koja ucitava parametre i jednacinu iz xml
int ucitajXML(vector<double>& koeficijenti, vector<double>& pocetni_uslovi, double& x_max) {
    cout << "\nParsiranje parametara za Runge Kutta metod: \n\n";

    xml_document doc;

    // load the XML file
    if (!doc.load_file("sample.xml")) return -1;

    xml_node tools = doc.child("RungeKutta").child("Parametri");

    //Za citanje pocetnih uslova
    for (xml_node_iterator it = tools.begin(); it != tools.end(); ++it) {
        cout << "Parametri: ";

        for (xml_attribute_iterator ait = it->attributes_begin(); ait != it->attributes_end(); ++ait) {
            pocetni_uslovi.push_back(std::stod(ait->value()));
            cout << " " << ait->name() << "=" << ait->value();
        }

        cout << endl;
    }

    xml_node tools_1 = doc.child("RungeKutta").child("Funkcija");
    std::string diferencijalna;


    //Za citanje diferencijalne jednacine
    for (xml_node_iterator it = tools_1.begin(); it != tools_1.end(); ++it) {
        cout << "Funkcija: ";

        for (xml_attribute_iterator ait = it->attributes_begin(); ait != it->attributes_end(); ++ait) {
            cout << " " << ait->name() << "=" << ait->value();
            diferencijalna = ait->value();
        }

        cout << endl;
    }
    koeficijenti = PronadjiKoeficijente(diferencijalna);

    xml_node tools_2 = doc.child("RungeKutta").child("Tacka_za_racunanje");


    //Za citanje tacke u kojoj se racuna
    for (xml_node_iterator it = tools_2.begin(); it != tools_2.end(); ++it) {
        cout << "x_max: ";

        for (xml_attribute_iterator ait = it->attributes_begin(); ait != it->attributes_end(); ++ait) {
            cout << " " << ait->name() << "=" << ait->value();
            x_max = std::stod(ait->value());
        }

        cout << endl;
    }

    cout << endl;

    // Ispisivanje dobivenih vrijednosti 
    cout << "Pocetni uslovi su: (x_0, y_0";
    for (int i = 2; i < pocetni_uslovi.size(); i++) {
        cout << ", y_0_" << i - 1;
        if (i == pocetni_uslovi.size() - 1) std::cout << ") = (";
    }
    for (int i = 0; i < pocetni_uslovi.size(); i++) {
        cout << pocetni_uslovi[i];
        if(i < pocetni_uslovi.size() - 1)
            cout << ", ";
    }
    cout << ")" << endl;
    cout << "Koeficijenti diferencijalne jednacine su (prvo uz x, zatim uz y): ";
    for (int i = 0; i < koeficijenti.size(); i++) {
        if (i == koeficijenti.size() - 1) std::cout << koeficijenti.at(i);
        else std::cout << koeficijenti.at(i) << ", ";
    }
    cout << endl;

    cout << "Tacka u kojoj se racuna funkcija je: x_max = " << x_max << endl;
    return 1;
}

//Funkcija koja vraca funkciju koja predstavlja diferencijalnu jednacinu na osnovu koeficijenata
function <double(double, double, vector<double>)> f(vector<double> koeficijenti) {
    return [koeficijenti](double x, double y, vector<double> u) {
        double result = koeficijenti[0] + koeficijenti[1] * x + koeficijenti[2] * y;
        for (int i = 0; i < u.size(); i++) {
            result += koeficijenti[i + 3] * u[i];
        }
        return result;
    };
}

//Funkcija koja izvrsava jedan korak RungeKutta4 za rjesavanje diferencijalnih jednacina prvog reda
double RK4Step(double f(double, double), double x, double y, double h) {
    std::vector<double> K(4);
    K.at(0) = f(x, y);
    K.at(1) = f(x + h / 2, y + h * K.at(0) / 2);
    K.at(2) = f(x + h / 2, y + h * K.at(1) / 2);
    K.at(3) = f(x + h / 2, y + h * K.at(2));
    return y + h * (K.at(0) + 2 * K.at(1) + 2 * K.at(2) + K.at(3)) / 6;
}

//RungeKutta 4 algoritam sa adaptacijom koraka za rjesavanje diferencijalnih jednacina prvog reda
double AdaptacijaKoraka(double f(double, double), double x_0, double y_0, double x_max, double h) {
    double x = x_0, y = y_0;
    double eps = 1e-5;
    while (x <= x_max) {
        double u = RK4Step(f, x, y, h / 2);
        double v = RK4Step(f, x + h / 2, u, h / 2);
        double w = RK4Step(f, x, y, h);
        double delta = std::fabs(w - v) / std::fabs(h);
        if (delta <= eps) {
            x = x + h;
            y = v;
        }
        h = h * fmin( 5, std::pow(0.9 * (eps / delta), 0.25) );
    }
    return y;
}

 //RungeKutta4 algoritam za rjesavanje diferencijalnih jednacina viseg reda
double RungeKutta(function <double(double, double, vector<double>)> f, double x_0, std::vector<double> y_0, double x_max, double h) {
    double x = x_0, y = y_0[0];
    vector<double> u, xx, yy;
    xx.push_back(x); yy.push_back(y);
    for (int i = 1; i < y_0.size(); i++) {
        u.push_back(y_0[i]);
    }

    vector<vector<double>> K(4, vector<double>(y_0.size())); 

    while (x <= x_max) {
        for (int k = 0; k < y_0.size() - 1; k++) {
            K[0][k] = u[k];
        }
        K[0][y_0.size() - 1] = f(x, y, u);
        for (int k = 0; k < y_0.size() - 1; k++) {
            K[1][k] = u[k] + K[0][k + 1] * h / 2.;
        }
        K[1][y_0.size() - 1] = f(x + h / 2, y + h * K[0][0] / 2, K[1]);
        for (int k = 0; k < y_0.size() - 1; k++) {
            K[2][k] = u[k] + K[1][k + 1] * h / 2.;
        }
        K[2][y_0.size() - 1] = f(x + h / 2, y + h * K[1][0] / 2, K[2]);
        for (int k = 0; k < y_0.size() - 1; k++) {
            K[3][k] = u[k] + K[2][k + 1] * h / 2.;
        }
        K[3][y_0.size() - 1] = f(x + h, y + h * K[2][0], K[3]);
        y += h * (K[0][0] + 2 * K[1][0] + 2 * K[2][0] + K[3][0]) / 6;
        for (int k = 0; k < y_0.size() - 1; k++) {
            u[k] += h * (K[0][k + 1] + 2 * K[1][k + 1] + 2 * K[2][k + 1] + K[3][k + 1]) / 6;
        }
        x += h;

        //vizualizacija

        xx.push_back(x); yy.push_back(y);
        plot.update(xx, yy);
        plt::pause(0.01);
    }
    return y;
}

int main(){
    vector<double> koeficijenti, pocetni_uslovi;
    double x_max;
    ucitajXML(koeficijenti, pocetni_uslovi, x_max);
    double x0 = pocetni_uslovi[0];
    pocetni_uslovi.erase(pocetni_uslovi.begin(), pocetni_uslovi.begin() + 1);
    std::vector<double> y0 = move(pocetni_uslovi);
    std::cout << "y0.size() = " << y0.size() << std::endl;
    plt::named_plot("okvir", std::vector<double>{0, 0, x_max, x_max, 0}, std::vector<double>{-10, 10, 10, -10, -10}, "w-");
    std::cout << RungeKutta(f(koeficijenti), x0, y0, x_max, 0.001);
    plt::show();
    return 0;
}
