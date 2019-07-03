//
// Created by kamilla on 03.04.19.
//

#ifndef FRESS_CPP_FRESS_H
#define FRESS_CPP_FRESS_H


#include <list>
#include <valarray>
#include <algorithm>
#include <random>
#include <tuple>
#include <vector>
#include <map>
#include <math.h>

typedef std::vector<int> Sequence_t, Sequence;
typedef std::pair<int, int> coord_t;
typedef std::vector<coord_t> Conformation_t, Conformation;
typedef std::map<int, std::vector<coord_t>>  contact_map_t ;
typedef std::map<coord_t, int>  Map_coordinate_to_int ;

class Lattice{

public:


    int ndim() { return 2; }
    int ndim2() {return 4;}
    int lattice_side;
    //contact_map_t map_of_contacts;
    //std:: map<int, std::pair<int, int>> map_int_to_coordinate;
    //Map_coordinate_to_int map_coordinate_to_int; //словарь для хранения всех соседей

    Lattice();
    Lattice(int seq_size );
    void create_lattice(int seq_size);
    int distance_lattice(coord_t point1, coord_t point2);

    inline std::vector<coord_t> get_contacts(coord_t point){
        return map_of_contacts[map_coordinate_to_int[point]];
    }

private:
    contact_map_t map_of_contacts;
    std:: map<int, std::pair<int, int>> map_int_to_coordinate;
    Map_coordinate_to_int map_coordinate_to_int; //словарь для хранения всех соседей

};

class Protein{
public:
    Sequence sequence;
   // int number_of_iterations;
    int E; // энергия текущей конформации
    int min_E;
    float T; //температура
    bool change_T; //чтобы корректно менять температуру каждый 50 000 шагов
    long int iterator; //для подсчета числа шагов
    std::vector<double> probabilities={}; //вектор вероятностей для выбора длины, вероятность пропорциональна 1/l
    Conformation conformation;
    std::vector<int> conformation_int;
    std::vector<Conformation> results;
    Lattice lattice;

    Protein();
    Protein(std::vector <int> sequence_input  );

    void find_minimum();
    void calculate_probabilities_for_l(int lmin = 2, int lmax = 12);
    int count_contacts();
    int count_contacts_dissected_t(Sequence_t &sequence, Conformation_t &conformation,  int t, int current_energy);
    int dissected(Sequence_t &sequence, Conformation_t &conformation);
    void regrowth_middle(int l, int start_position);
    void regrowth_start(int l );
    void regrowth_end(int l );


private:




};









#endif //FRESS_CPP_FRESS_H
