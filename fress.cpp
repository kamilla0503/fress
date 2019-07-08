//
// Created by kamilla on 03.04.19.
//

#include "fress.h"
#include <iostream>
#include <fstream>

Lattice::Lattice() {};

Lattice::Lattice(int seq_size ) {
    lattice_side = seq_size+3;
    int l_s = seq_size+2; //чтобы было точно достатчно места при граничных условия для тора
    //создается словарь, в котором ключ - номер точки на квардатной решетке
    // значения - номера четырех соседних точек
    int x, y;
    div_t n;
    for (int i =0; i<(l_s+1)*(l_s+1); i++){
        //map_int_to_coordinate[i]=std::make_pair(i%(l_s+1), i /(l_s+1)  );
        //map_coordinate_to_int[map_int_to_coordinate[i]] = i;
        map_of_contacts_int[i] = {};
        map_of_contacts_int[i].push_back(i+1);
        map_of_contacts_int[i].push_back(i-1);
        map_of_contacts_int[i].push_back(i+l_s+1);
        map_of_contacts_int[i].push_back(i-l_s-1);
        n=div(i, l_s+1);
        x=n.rem;
        y=n.quot;
        for (int j =0; j<ndim2(); j++){
            if(x==0){
                map_of_contacts_int[i][1] = i+l_s;
            }
            if(x==l_s){
                map_of_contacts_int[i][0] = i-l_s;
            }
            if(y==0){
                map_of_contacts_int[i][3] = l_s*(l_s+1)+x;
            }
            if(y==l_s){
                map_of_contacts_int[i][2] = x;
            }

        }

    }

};

void Lattice::create_lattice(int seq_size){
    lattice_side = seq_size+3;
    int l_s = seq_size+2; //чтобы было точно достатчно места при граничных условия для тора
    //создается словарь, в котором ключ - номер точки на квардатной решетке
    // значения - номера четырех соседних точек
    int x, y;
    div_t n;
    for (int i =0; i<(l_s+1)*(l_s+1); i++){
        //map_int_to_coordinate[i]=std::make_pair(i%(l_s+1), i /(l_s+1)  );
        //map_coordinate_to_int[map_int_to_coordinate[i]] = i;
        map_of_contacts_int[i] = {};
        map_of_contacts_int[i].push_back(i+1);
        map_of_contacts_int[i].push_back(i-1);
        map_of_contacts_int[i].push_back(i+l_s+1);
        map_of_contacts_int[i].push_back(i-l_s-1);
        n=div(i, l_s+1);
        x=n.rem;
        y=n.quot;
        for (int j =0; j<ndim2(); j++){
            if(x==0){
                map_of_contacts_int[i][1] = i+l_s;
            }
            if(x==l_s){
                map_of_contacts_int[i][0] = i-l_s;
            }
            if(y==0){
                map_of_contacts_int[i][3] = l_s*(l_s+1)+x;
            }
            if(y==l_s){
                map_of_contacts_int[i][2] = x;
            }

        }

    }





}

Protein::Protein(){};

Protein::Protein(std::vector<int> sequence_input ) {
    sequence=sequence_input;
    iterator=0;
    T=3.5;
    change_T=false;
    calculate_probabilities_for_l();
   lattice.create_lattice(sequence.size());
    for (int i=0; i<sequence.size(); i++){
        conformation.push_back(i);
    }
    E = count_contacts();
    min_E = E; //сохраняем минимальное найденное значение
    results={};
}

void Protein::calculate_probabilities_for_l(int lmin, int lmax) {
    //создается вектор с вероятностями, с которыми будут выбираться l
    double one_over_l = 0.0;
    for (int i=lmin; i <=lmax; i++ ){
        one_over_l=one_over_l+1.0/i;
    }
    double  k = 1/one_over_l;
    for (int i=lmin; i <=lmax; i++ ){
        probabilities.push_back(k/i);
    }
    for (int i =1; i<probabilities.size(); i++){
        probabilities[i]=probabilities[i]+probabilities[i-1];
    }
}

int Protein::count_contacts(){
    //считается число топологических контактов НН
    int hh = 0;
    int position;
    for (int i =1; i<sequence.size()-1; i++){
        for ( coord_t step : lattice.get_contacts(conformation[i]) ){
            if ( step!=conformation[i-1] && step!=conformation[i+1] && std::find(conformation.begin(), conformation.end(), step) !=conformation.end()  ){
                position=std::distance(conformation.begin(),find(conformation.begin(), conformation.end(), step));
                hh=hh+sequence[i]*sequence[position];
            }
        }
    }
    //концы конформации обрабатываются отдельно
    for ( coord_t step : lattice.get_contacts(conformation[0]) ) {
        if (step != conformation[1] &&
            std::find(conformation.begin(), conformation.end(), step) != conformation.end()) {
            position = std::distance(conformation.begin(), find(conformation.begin(), conformation.end(), step));
            hh = hh + sequence.front() * sequence[position];
        }
    }
    for ( coord_t step : lattice.get_contacts(conformation.back()) ){
        if ( step!=conformation[conformation.size()-2]  && std::find(conformation.begin(), conformation.end(), step) !=conformation.end()  ){
            position=std::distance(conformation.begin(),find(conformation.begin(), conformation.end(), step));
            hh=hh+sequence.back()*sequence[position];
        }
    }
    return  (-1*div(hh, 2).quot);
}


int Lattice::distance_lattice(coord_t point_1, coord_t point_2){
    //функция для подсчета расстояния между двумя точками на плоскости с учетом граничных условий
    div_t point1, point2;
    point1 = div(point_1, lattice_side);
    point2 = div(point_2, lattice_side);
    int p1 = abs(point2.rem-point1.rem);
    int p2 = lattice_side -  abs(point2.rem-point1.rem);
    int part1 = std::min(p1, p2);
    p1 = abs(point2.quot-point1.quot);
    p2 = lattice_side -  abs(point2.quot-point1.quot  );
    int part2 = std::min(p1, p2);
    return  part1 + part2;
}

void Protein::regrowth_end(int l ){
    //для выращивания участка длины l на конце цепочки
    Conformation_t  C_t ;
    Sequence_t seq_t;
    C_t.reserve(sequence.size());
    C_t.resize(sequence.size()-l);
    std::copy(conformation.begin(), conformation.begin()+sequence.size()-l, C_t.begin());
    seq_t.reserve(sequence.size());
    seq_t.resize(sequence.size()-l);
    std::copy(sequence.begin(), sequence.begin()+sequence.size()-l, seq_t.begin());
    int current_energy = dissected(seq_t, C_t );
    int temp_e;
    static std::vector<std::pair<float, int>> probabilities_to_move;
    probabilities_to_move.resize(4, std::make_pair(0.0, 0));
    std::vector<int> energies;
    energies.resize(lattice.ndim2(), 0);
    float sum_probabilities=0.0;
    seq_t.emplace_back( sequence[sequence.size()-l]);
    for (int i =0; i<lattice.ndim2(); i++ ){
        //проверка на самопересечения + проверка первого шага, чтобы не начать строить в том же направлении, что и в старой конформации
        if (lattice.get_contacts(C_t.back())[i] == conformation[sequence.size()-l] ||  std::find(C_t.begin(), C_t.end(), lattice.get_contacts(C_t.back())[i]   )!=C_t.end()  ){
            //так как точка не подходит, то вероятность попасть в нее равна нулю
            probabilities_to_move[i] = std::make_pair(0.0, i);
            energies[i] = 0;//strange, but it for time economy
            continue;
        }
        else{
            //точка подходит
            //текущая  функция энергия должна работать корректно
            C_t.emplace_back(lattice.get_contacts(C_t.back())[i]);
            temp_e = count_contacts_dissected_t(seq_t, C_t,  sequence.size() - l,
                                                current_energy);
            energies[i] = temp_e;
            probabilities_to_move[i] = std::make_pair( exp(-(temp_e-current_energy)/T), i);
            // это все  energy-guided
            C_t.pop_back();
        }
        sum_probabilities=sum_probabilities+probabilities_to_move[i].first;
    }
    if(sum_probabilities==0.0){
     // значит двигаться некуда, попытка перестройки оказалась неуспешной
        return;
    }
    sort(probabilities_to_move.begin(), probabilities_to_move.end());
    for (int i =0; i<lattice.ndim2(); i++ ){
        probabilities_to_move[i].first= probabilities_to_move[i].first/sum_probabilities;
        // чтобы вероятности стали вероятностями
    }
    for (int i =1; i<lattice.ndim2(); i++ ){
        probabilities_to_move[i].first= probabilities_to_move[i].first + probabilities_to_move[i-1].first;
        // чтобы осуществить выбор с нужными вероятностями
    }
    std::default_random_engine generator(std::random_device{}() );
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    double q = distribution(generator); // случайно число от 0 до 1
    for (int i =0; i<lattice.ndim2(); i++ ){
        //забираем соответствующую координату
        if (q<probabilities_to_move[i].first){
            C_t.emplace_back( lattice.get_contacts(C_t.back())[probabilities_to_move[i].second]  );
            current_energy = energies[probabilities_to_move[i].second];
            break;
        }
    }
    // выстраиваем остальную часть цепочки
    for (int t = sequence.size()-l+1; t<sequence.size(); t++) {
        seq_t.emplace_back( sequence[t]);
        sum_probabilities = 0.0;
        for (int i =0; i<lattice.ndim2(); i++ ){
            // тут проверка только на самопересечения
            if(std::find(C_t.begin(), C_t.end(),lattice.get_contacts(C_t.back())[i])==C_t.end() ){
                C_t.emplace_back(lattice.get_contacts(C_t.back())[i]);
                temp_e = count_contacts_dissected_t(seq_t, C_t, t,
                                                    current_energy);
                energies[i] = temp_e;
                probabilities_to_move[i] = std::make_pair( exp(-(temp_e-current_energy)/T), i);
                C_t.pop_back();
            }
            else {
                //точка не подходит
                probabilities_to_move[i] = std::make_pair(0.0, i);
                energies[i] = 0;
            }
            sum_probabilities=sum_probabilities+probabilities_to_move[i].first;
        }
        if(sum_probabilities==0.0){
            return;
        }
        sort(probabilities_to_move.begin(), probabilities_to_move.end());
        for (int i =0; i<lattice.ndim2(); i++ ){
            probabilities_to_move[i].first= probabilities_to_move[i].first/sum_probabilities;
        }
        for (int i =1; i<lattice.ndim2(); i++ ){
            probabilities_to_move[i].first= probabilities_to_move[i].first + probabilities_to_move[i-1].first;
        }
        q = distribution(generator);
        for (int i =0; i<lattice.ndim2(); i++ ){
            if (q<probabilities_to_move[i].first){
                C_t.emplace_back(lattice.get_contacts(C_t.back())[probabilities_to_move[i].second]  );
                current_energy = energies[probabilities_to_move[i].second];
                break;
            }
        }
    }
// дополнительная проверка, что раз мы дошли до этой части, участок нужной длины достроен
// по сути не нужно проверять, но мне так спокойнее
    if( C_t.size()==sequence.size()  ){
        current_energy = dissected(sequence, C_t); // пересчет
        q =distribution(generator);
        //это как раз вероятность, с которой будет приниматься новая конфромация
        // в оригинале это p=min( 1,   exp(..)  )
        // но так как q - число от 0 до 1, можно обойтись без функции min
        float probability_to_accept = exp(-(current_energy-E)/T);
        if(q<probability_to_accept){
            E = current_energy;
            conformation=C_t;
            if(E<min_E){
                results.clear();
                results.emplace_back(conformation);
                min_E=E;
                std:: cout << "New minimum fonded: " << min_E << std::endl ;
                iterator=iterator+1;
                if(change_T== false){
                    change_T=true;
                }
                std::ofstream out;
                out.open("coordinates_first_result.txt");
                for (int c : conformation) {
                  //  out << c.first << " " << c.second << "   ";
                    out << c;
                    out << std::endl;
                }
                out.close();
            }
            if (E==min_E){
                results.emplace_back(conformation);
            }
        }
    }
    else{
        //сюда не должны попадать никогда
        std::cout << " fail int the end " << std:: endl;
        std:: cout << "new size " << C_t.size() << std:: endl;
        return ;
    }
}

void Protein::regrowth_start(int l ) {
    Conformation_t  C_t ;
    Sequence_t seq_t;
    C_t.reserve(sequence.size());
    C_t.resize(sequence.size()-l);
    std::copy(conformation.begin()+l, conformation.end(), C_t.begin());
    seq_t.reserve(sequence.size());
    seq_t.resize(sequence.size()-l);
    std::copy(sequence.begin()+l, sequence.end(), seq_t.begin());
    int current_energy = dissected(seq_t, C_t);
    int temp_e;
    static std::vector <std::pair<float, int>> probabilities_to_move;
    probabilities_to_move.resize(4, std::make_pair(0.0, 0));
    std::vector<int> energies;
    energies.resize(lattice.ndim2(), 0);
    float sum_probabilities=0.0;
    seq_t.insert(seq_t.begin(), sequence[l-1]);
    for (int i =0; i<lattice.ndim2(); i++ ){
        //убираю вариант старой точки и самопересечения
        if (lattice.get_contacts(C_t[0])[i] == conformation[l-1] ||  std::find(C_t.begin(), C_t.end(), lattice.get_contacts(C_t[0])[i]   )!=C_t.end()  ){
            probabilities_to_move[i] = std::make_pair(0.0, i);
            energies[i] = 0;//strange, but it for time economy
            continue;
        }
        else{
            //точка подходит
            // Лучше потом переделать функцию для энергии
            C_t.insert(C_t.begin(), lattice.get_contacts(C_t[0])[i]  );
            temp_e = count_contacts_dissected_t(seq_t, C_t,   0, current_energy);
            energies[i] = temp_e;
            probabilities_to_move[i] = std::make_pair( exp(-(temp_e-current_energy)/T), i);
            C_t.erase(C_t.begin());
        }
        sum_probabilities=sum_probabilities+probabilities_to_move[i].first;
    }
    if(sum_probabilities==0.0){
        //неудачная попытка
        return;
    }
    sort(probabilities_to_move.begin(), probabilities_to_move.end());
    for (int i =0; i<lattice.ndim2(); i++ ){
        probabilities_to_move[i].first= probabilities_to_move[i].first/sum_probabilities;
    }
    for (int i =1; i<lattice.ndim2(); i++ ){
        probabilities_to_move[i].first= probabilities_to_move[i].first + probabilities_to_move[i-1].first;
    }
    std::default_random_engine generator(std::random_device{}() );
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    double q = distribution(generator);
    for (int i =0; i<lattice.ndim2(); i++ ){
        if (q<probabilities_to_move[i].first){
            C_t.insert(C_t.begin(),lattice.get_contacts(C_t[0])[probabilities_to_move[i].second]  );
            current_energy = energies[probabilities_to_move[i].second];
            break;
        }
    }
    for (int t = l-2; t>-1; t--){
        seq_t.insert(seq_t.begin(), sequence[t]);
        sum_probabilities = 0.0;
        for (int i =0; i<lattice.ndim2(); i++ ){
            if(std::find(C_t.begin(), C_t.end(), lattice.get_contacts(C_t[0])[i])==C_t.end() ){
                C_t.insert(C_t.begin(),lattice.get_contacts(C_t[0])[i]  );
                temp_e = count_contacts_dissected_t(seq_t, C_t, 0,
                                                    current_energy);
                energies[i] = temp_e;
                probabilities_to_move[i] = std::make_pair( exp(-(temp_e-current_energy)/T), i);
                C_t.erase(C_t.begin());
            }
            else {
                //точка не подходит
                probabilities_to_move[i] = std::make_pair(0.0, i);
                energies[i] = 0;
            }
            sum_probabilities=sum_probabilities+probabilities_to_move[i].first;
        }
        if(sum_probabilities==0.0){
            return;
        }
        sort(probabilities_to_move.begin(), probabilities_to_move.end());
        for (int i =0; i<lattice.ndim2(); i++ ){
            probabilities_to_move[i].first= probabilities_to_move[i].first/sum_probabilities;
        }
        for (int i =1; i<lattice.ndim2(); i++ ){
            probabilities_to_move[i].first= probabilities_to_move[i].first + probabilities_to_move[i-1].first;
        }
        q = distribution(generator);
        for (int i =0; i<lattice.ndim2(); i++ ){
            if (q<probabilities_to_move[i].first){
                C_t.insert(C_t.begin(),lattice.get_contacts(C_t[0])[probabilities_to_move[i].second]  );
                current_energy = energies[probabilities_to_move[i].second];
                break;
            }
        }
    }
    if( C_t.size()==sequence.size()  ){
        current_energy = dissected(sequence, C_t);
        q =distribution(generator);
        float probability_to_accept = exp(-(current_energy-E)/T);
        if(q<probability_to_accept){
            E = current_energy;
            conformation=C_t;
            iterator=iterator+1;
            if(change_T== false){
                change_T=true;
            }
            if(E<min_E){
                results.clear();
                results.emplace_back(conformation);
                min_E=E;
                std:: cout << "New minimum founded: " << min_E << std::endl ;
                std::ofstream out;
                out.open("coordinates_first_result.txt");
                for (int c : conformation) {
                    //out << c.first << " " << c.second << "   ";
                    out << c ;
                    out << std::endl;
                }
                out.close();
            }
            if (E==min_E){
                results.emplace_back(conformation);
            }
        }
    }
    else{
        std::cout << " fail int the end " << std:: endl;
        std:: cout << "new size " << C_t.size() << std:: endl;
        return ;
    }
}


void Protein::regrowth_middle(int l, int start_position){
    int end_position = start_position+l-1;
    Conformation_t C_t, C_t_temp ;
    Sequence_t seq_t, seq_t_temp;
    //std::copy(sequence.begin(), seq_t.begin()+ start_position, C_t_temp.begin());
    C_t.reserve(sequence.size());
    C_t_temp.resize(start_position);
    std::copy(conformation.begin(), conformation.begin()+ start_position, C_t_temp.begin());
    C_t = C_t_temp;
    C_t_temp.resize(sequence.size()-start_position-l);
    std::copy(conformation.begin()+end_position+1, conformation.end(), C_t_temp.begin());
    C_t.insert(C_t.end(), C_t_temp.begin(), C_t_temp.end());
    C_t_temp.clear();
    //это удаление куска конформации с start_position до end_position
    seq_t.reserve(sequence.size());
    seq_t_temp.resize(start_position);
    std::copy(sequence.begin(), sequence.begin()+ start_position, seq_t_temp.begin());
    seq_t = seq_t_temp;
    seq_t_temp.resize(sequence.size()-start_position-l);
    //seq_t_temp.clear();
    std::copy(sequence.begin()+end_position+1, sequence.end(), seq_t_temp.begin());
   // seq_t.reserve(seq_t.size()+seq_t_temp.size() );
    seq_t.insert(seq_t.end(), seq_t_temp.begin(), seq_t_temp.end());
    seq_t_temp.clear();
    //это удаление куска последовательности
    int current_energy = dissected(seq_t, C_t );
    seq_t.insert(seq_t.begin()+start_position, sequence[start_position]);
    int temp_e;
    static std::vector <std::pair<float, int>> probabilities_to_move;
    probabilities_to_move.resize(4, std::make_pair(0.0, 0));
    std::vector<int> energies;
    energies.resize(lattice.ndim2(), 0);
    std::pair<int, int> point;
    float sum_probabilities=0.0;
    for (int i =0; i<lattice.ndim2(); i++ ){
         if (lattice.get_contacts(C_t[start_position-1])[i] == conformation[start_position]){
            probabilities_to_move[i] = std::make_pair(0.0, i);
            energies[i] = 0;//strange, but it is for time economy
            continue;
        }
        else  if( std::find(C_t.begin(), C_t.end(), lattice.get_contacts(C_t[start_position-1])[i]   )==C_t.end() && lattice.distance_lattice(lattice.get_contacts(C_t[start_position-1])[i], conformation[end_position+1]) <=abs(end_position+1-start_position)  ){
            // Лучше потом переделать функцию для энергии
            C_t.insert(C_t.begin()+start_position,lattice.get_contacts(C_t[start_position-1])[i]  );

            energies[i] = count_contacts_dissected_t(seq_t, C_t,  start_position,
                                                     current_energy);
             probabilities_to_move[i] = std::make_pair( exp(-(temp_e-current_energy)/T), i);
            C_t.erase(C_t.begin()+start_position);
        }
        else{
            probabilities_to_move[i] = std::make_pair(0.0, i);
            energies[i] = 0;//strange, but it for time economy
            continue;
        }
        sum_probabilities=sum_probabilities+probabilities_to_move[i].first;
    }
    if(sum_probabilities==0.0){
        return;
    }
    sort(probabilities_to_move.begin(), probabilities_to_move.end());
    for (int i =0; i<lattice.ndim2(); i++ ){
        probabilities_to_move[i].first= probabilities_to_move[i].first/sum_probabilities;
    }
    for (int i =1; i<4; i++ ){
        probabilities_to_move[i].first= probabilities_to_move[i].first + probabilities_to_move[i-1].first;
    }
    std::default_random_engine generator(std::random_device{}() );
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    double q = distribution(generator);
    for (int i =0; i<lattice.ndim2(); i++ ){
        if (q<probabilities_to_move[i].first){
            C_t.insert(C_t.begin()+start_position,lattice.get_contacts(C_t[start_position-1])[probabilities_to_move[i].second]  );
            current_energy = energies[probabilities_to_move[i].second];
            break;
        }
    }
    for (int t = start_position+1; t< end_position+1; t++){
        seq_t.insert(seq_t.begin()+t, sequence[t]);
        sum_probabilities = 0.0;
        for (int i =0; i<lattice.ndim2(); i++ ){
            auto tq=lattice.get_contacts(C_t[t-1]) ;
            if(std::find(C_t.begin(), C_t.end(), lattice.get_contacts(C_t[t-1])[i])==C_t.end() &&  lattice.distance_lattice( lattice.get_contacts(C_t[t-1])[i] ,conformation[end_position+1] )<=abs(end_position+1-t)){
                C_t.insert(C_t.begin()+t,lattice.get_contacts(C_t[t-1])[i]  );
                temp_e = count_contacts_dissected_t(seq_t, C_t,  t,
                                                    current_energy);
                energies[i] = temp_e;
                probabilities_to_move[i] = std::make_pair( exp(-(temp_e-current_energy)/T), i);
                C_t.erase(C_t.begin()+t);
            }
            else {
                //точка подходит
                probabilities_to_move[i] = std::make_pair(0.0, i);
                energies[i] = 0;
            }
            sum_probabilities=sum_probabilities+probabilities_to_move[i].first;
        }
        if(sum_probabilities==0.0){
            return;
        }
        sort(probabilities_to_move.begin(), probabilities_to_move.end());
        for (int i =0; i<lattice.ndim2(); i++ ){
            probabilities_to_move[i].first= probabilities_to_move[i].first/sum_probabilities;
        }
        for (int i =1; i<lattice.ndim2(); i++ ){
            probabilities_to_move[i].first= probabilities_to_move[i].first + probabilities_to_move[i-1].first;
        }
        q = distribution(generator);
        for (int i =0; i<lattice.ndim2(); i++ ){
            if (q<probabilities_to_move[i].first){
                C_t.insert(C_t.begin()+t,lattice.get_contacts(C_t[t-1])[probabilities_to_move[i].second]  );
                current_energy = energies[probabilities_to_move[i].second];
                break;
            }
        }
    }
// проверка по сути не нужна, но мне так спокойнее
    if( C_t.size()==sequence.size()  ){
        current_energy = dissected(sequence, C_t );
        q =distribution(generator);
        float probability_to_accept = exp(-(current_energy-E)/T);
        if(q<probability_to_accept){
            E = current_energy;
            conformation=C_t;
            if(E<min_E){
                results.clear();
                results.emplace_back(conformation);
                min_E=E;
                iterator=iterator+1;
                if(change_T== false){
                    change_T=true;
                }
                std:: cout << "New minimum fonded: " << min_E << std::endl ;
                std::ofstream out;
                out.open("coordinates_first_result.txt");
                for (auto c : conformation) {
                  //  out << c.first << " " << c.second << "   ";
                  out << c;
                    out << std::endl;
                }
                out.close();
            }
            if (E==min_E){
                results.emplace_back(conformation);
            }
        }
    }
    else{
    std::cout << " fail int the end " << std:: endl;
    std:: cout << "new size " << C_t.size() << std:: endl;
    return ;
    }
}

void Protein::find_minimum() {
    double q;
    int l, start_position;
    std::default_random_engine generator(std::random_device{}() );
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    std::default_random_engine generators1(std::random_device{}() );
    for (int iteration =1; iteration<1000000000000; iteration++) {
        q = distribution(generator);
        //std::cout<< q << " " ;
        for (int i = 0; i < probabilities.size(); i++) {
            if (q < probabilities[i]) {

                l =i+2;
                break;
            }
        }
        std::uniform_int_distribution<int> distribution1(0,sequence.size()-l);
        start_position = distribution1(generators1);
        if(start_position==0){
            regrowth_start(l);
        }
        else if (start_position==sequence.size()-l){
            regrowth_end(l);
        }
        else{
            regrowth_middle(l, start_position);
        }
        //std:: cout << l << " " ;
        if( iterator%50000==0 && change_T && T>0.1){
            T=T*0.98;
            std:: cout << "50000 moves are made; new T : " << T << std:: endl;
            change_T = false;
        }
        if(iteration%10000==0){
            std::cout << "Number of attempts :  " << iteration<< std::endl;
        }
        if(min_E==-14){
            std::cout << "-9 achieved! " << std:: endl;
            break;
        }
    }
    std::cout << std::endl;
};


int Protein::count_contacts_dissected_t(Sequence_t &sequence1, Conformation_t &conformation1, int t,
                                        int current_energy) {
    int new_energy = current_energy;
    int position;
    if(t!=0 && t!=sequence.size()-1){
        for (coord_t step : lattice.get_contacts(conformation1[t])) {
            if (step != conformation1[t - 1] && step != conformation1[t + 1] &&
                std::find(conformation1.begin(), conformation1.end(), step) != conformation1.end()) {
                position = std::distance(conformation1.begin(), find(conformation1.begin(), conformation1.end(), step));
                new_energy = new_energy - sequence1[t] * sequence1[position];
            }
        }
    }
    else if(t==0){
        for (coord_t step : lattice.get_contacts(conformation1[0])) {
            if (  step != conformation1[  1] &&
                std::find(conformation1.begin(), conformation1.end(), step) != conformation1.end()) {
                position = std::distance(conformation1.begin(), find(conformation1.begin(), conformation1.end(), step));
                new_energy = new_energy - sequence1[0] * sequence1[position];
            }
        }
    }
    else{
        for (coord_t step : lattice.get_contacts(conformation1[t])) {
            if (step != conformation1[conformation1.size()-2]  &&
                std::find(conformation1.begin(), conformation1.end(), step) != conformation1.end()) {
                position = std::distance(conformation1.begin(), find(conformation1.begin(), conformation1.end(), step));
                new_energy = new_energy - sequence1[t] * sequence1[position];
            }
        }
    }
    return new_energy ;
}

int Protein::dissected(Sequence_t &sequence1, Conformation_t &conformation1) {
    int hh = 0;
    int position;
    for (int i = 1; i < sequence1.size()-1; i++) {
        for (coord_t step : lattice.get_contacts(conformation1[i])) {
            if (step != conformation1[i - 1] && step != conformation1[i + 1] &&
                std::find(conformation1.begin(), conformation1.end(), step) != conformation1.end()) {
                position = std::distance(conformation1.begin(), find(conformation1.begin(), conformation1.end(), step));
                hh = hh + sequence1[i] * sequence1[position];
            }
        }
    }
    for (coord_t step : lattice.get_contacts(conformation1[0])) {
        if (step != conformation1[1] &&
            std::find(conformation1.begin(), conformation1.end(), step) != conformation1.end()) {
            position = std::distance(conformation1.begin(), find(conformation1.begin(), conformation1.end(), step));
            hh = hh + sequence1[0] * sequence1[position];
        }
    }
    for (coord_t step : lattice.get_contacts(conformation1.back())) {
        if (step != conformation1[conformation1.size() - 2] &&
            std::find(conformation1.begin(), conformation1.end(), step) != conformation1.end()) {
            position = std::distance(conformation1.begin(), find(conformation1.begin(), conformation1.end(), step));
            hh = hh + sequence1.back() * sequence1[position];
        }
    }
    return (-1 * div(hh, 2).quot);
}