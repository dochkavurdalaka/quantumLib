#pragma once
#pragma once
#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>


class Gate {
private:
    std::vector<std::vector<std::complex<double>>> data;
    size_t size;

public:
    Gate(std::vector<std::vector<std::complex<double>>> temp) :data(temp), size(temp.size()) {
    }


    std::vector<std::complex<double>>& operator[](int index) {
        return data[index];
    }

    std::vector<std::complex<double>> operator[](int index)  const {
        return data[index];
    }
};


Gate IGate = Gate({ {std::complex<double>(1,0),std::complex<double>(0,0)} , {std::complex<double>(0,0),std::complex<double>(1,0)} });
Gate XGate = Gate({ {std::complex<double>(0,0),std::complex<double>(1,0)} , {std::complex<double>(1,0),std::complex<double>(0,0)} });
Gate HGate = Gate({ {std::complex<double>(1.0 / sqrt(2.0),0),std::complex<double>(1.0 / sqrt(2.0),0)} , {std::complex<double>(1.0 / sqrt(2.0),0),std::complex<double>(-1.0 / sqrt(2.0),0)} });
Gate ZGate = Gate({ {std::complex<double>(1,0),std::complex<double>(0,0)} , {std::complex<double>(0,0),std::complex<double>(-1,0)} });


class State {
private:
    std::vector<std::complex<double>> data;
    size_t size;
    size_t num;


    static unsigned int insert_bit(unsigned int t, int bit) {
        int size = sizeof(t) * 8;
        unsigned int temp1;
        if (bit == 0) {
            temp1 = 0;
        }
        else {
            temp1 = (t << (size - bit)) >> (size - bit);
        }
        unsigned int temp2 = (t >> bit) << (bit + 1);
        return (temp1 | temp2);
    }



public:
    State(std::vector<std::complex<double>> temp) :data(temp), size(temp.size()), num((int)round(log2(temp.size()))) {

    }

    std::complex<double>& operator[](int index) {
        return data[index];
    }

    std::complex<double> operator[](int index)  const {
        return data[index];
    }

    size_t getNum() {
        return num;
    }

    void MeasureAll() {

        double alpha = (double)rand() / RAND_MAX;

        double S = 0.0;
        int i = 0;
        while (alpha > S) {
            S += (data[i] * std::conj(data[i])).real();
            i++;
        }

        for (int i = 0; i < size; i++)
            data[i] = 0;

        data[i - 1] = 1.0;
    }

    void MeasureBit(int bit) {
        std::vector<std::complex<double>> vector_temp(size, std::complex<double>(0, 0));
        double S0 = 0.0;
        double S1 = 0.0;

        int temp = 1u << bit;
        for (int i = 0; i < data.size(); i += temp * 2) {
            for (int j = i; j < i + temp; j++) {
                S0 += (data[j] * std::conj(data[j])).real();
            }
            for (int j = i + temp; j < i + 2 * temp; j++) {
                S1 += (data[j] * std::conj(data[j])).real();
            }
        }

        double alpha = (double)rand() / RAND_MAX;
        if (alpha <= S0) {
            for (int i = 0; i < data.size(); i += temp * 2) {
                for (int j = i; j < i + temp; j++) {
                    data[j] /= sqrt(S0);
                }
            }
        }
        else {
            for (int i = 0; i < data.size(); i += temp * 2) {
                for (int j = i + temp; j < i + 2 * temp; j++) {
                    data[j] /= sqrt(S1);
                }
            }
        }

    }


    void CX(int i, int j) {
        for (int k = 0; k < size / 4; k++) {
            unsigned int temper = (unsigned int)k;
            if (i < j) {
                temper = insert_bit(temper, i);
                temper = insert_bit(temper, j);
            }
            else {
                temper = insert_bit(temper, j);
                temper = insert_bit(temper, i);
            }

            temper = temper | (1u << i);

            std::complex<double> temp = data[temper];
            data[temper] = data[temper | (1u << j)];
            data[temper | (1u << j)] = temp;
        }
    }


    void CP(double fi, int i, int j ) {
        for (int k = 0; k < size / 4; k++) {
            unsigned int temper = (unsigned int)k;
            if (i < j) {
                temper = insert_bit(temper, i);
                temper = insert_bit(temper, j);
            }
            else {
                temper = insert_bit(temper, j);
                temper = insert_bit(temper, i);
            }

            temper = temper | (1u << i) | (1u << j);

            data[temper] = data[temper] * std::complex<double>(cos(fi), sin(fi));
        }
    }

    void CGate(Gate gate, int control, int target) {
        for (int k = 0; k < size / 4; k++) {
            unsigned int temper = (unsigned int)k;
            if (control < target) {
                temper = insert_bit(temper, control);
                temper = insert_bit(temper, target);
            }
            else {
                temper = insert_bit(temper, target);
                temper = insert_bit(temper, control);
            }

            unsigned int one = temper | (1u << control);
            unsigned int two = temper | (1u << control) | (1u << target);

            std::vector<std::complex<double>> temp(2, std::complex<double>(0, 0));

            temp[0] = data[one] * gate[0][0] + data[two] * gate[0][1];
            temp[1] = data[one] * gate[1][0] + data[two] * gate[1][1];

            data[one] = temp[0];
            data[two] = temp[1];
        }
    }
    
    void CCGate(Gate gate, int contr1, int contr2, int target) {
        for (int k = 0; k < size / 8; k++) {
            unsigned int temper = (unsigned int)k;
            std::vector <int> sort_temp = { contr1, contr2, target};
            std::sort(sort_temp.begin(), sort_temp.end());
            temper = insert_bit(temper, sort_temp[0]);
            temper = insert_bit(temper, sort_temp[1]);
            temper = insert_bit(temper, sort_temp[2]);


            unsigned int one = temper | (1u << contr1) | (1u << contr2);
            unsigned int two = temper | (1u << contr1) | (1u << contr2) | (1u << target);


            std::vector<std::complex<double>> temp(2, std::complex<double>(0, 0));

            temp[0] = data[one] * gate[0][0] + data[two] * gate[0][1];
            temp[1] = data[one] * gate[1][0] + data[two] * gate[1][1];

            data[one] = temp[0];
            data[two] = temp[1];
        }
    }
    

    void CCX(int i, int j, int o) {
        for (int k = 0; k < size / 8; k++) {
            unsigned int temper = (unsigned int)k;
            std::vector <int> sort_temp = { i, j, o };
            std::sort(sort_temp.begin(), sort_temp.end());
            temper = insert_bit(temper, sort_temp[0]);
            temper = insert_bit(temper, sort_temp[1]);
            temper = insert_bit(temper, sort_temp[2]);


            temper = temper | (1u << i);
            temper = temper | (1u << j);


            std::complex<double> temp = data[temper];
            data[temper] = data[temper | (1u << o)];
            data[temper | (1u << o)] = temp;
        }
    }

    void CCP(double fi, int i, int j, int o) {
        for (int k = 0; k < size / 8; k++) {
            unsigned int temper = (unsigned int)k;
            std::vector <int> sort_temp = { i, j, o };
            std::sort(sort_temp.begin(), sort_temp.end());
            temper = insert_bit(temper, sort_temp[0]);
            temper = insert_bit(temper, sort_temp[1]);
            temper = insert_bit(temper, sort_temp[2]);


            temper = temper | (1u << i) | (1u << j) | (1u << o);


            data[temper] = data[temper] * std::complex<double>(cos(fi), sin(fi));
        }
    }



    void apply(Gate gate, int bit) {
        std::vector<std::complex<double>> t(2, std::complex<double>(0, 0));
        int two = 1u << bit;
        for (int i = 0; i < size; i += two * 2) {
            for (int j = i; j < i + two; j++) {
                t[0] = data[j] * gate[0][0] + data[j + two] * gate[0][1];
                t[1] = data[j] * gate[1][0] + data[j + two] * gate[1][1];

                data[j] = t[0];
                data[j + two] = t[1];
            }
        }
    }

    void I(int bit) {
        apply(IGate, bit);
    }

    void X(int bit) {
        apply(XGate, bit);
    }

    void H(int bit) {
        apply(HGate, bit);
    }

    void Z(int bit) {
        apply(ZGate, bit);
    }

    void P(int bit, double fi) {
        Gate P_ = Gate({ {std::complex<double>(1,0),std::complex<double>(0,0)} , {std::complex<double>(0,0),std::complex<double>(cos(fi),sin(fi))} });
        apply(P_, bit);
    }



    friend std::ostream& operator<< (std::ostream& out, const State& other);

};





std::ostream& operator<< (std::ostream& out, const State& other)
{


    for (size_t i = 0; i < other.size; i++) {
        out << other[i] << "\n";
    }


    return out;
}
