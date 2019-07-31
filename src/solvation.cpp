/****************************************************************************
 *   Copyright (C) 2017   Romain Deniau                                     *
 *   Copyright (C) 2017   Jean-Christophe Gelly                             *
 *   Copyright (C) 2017   Yassine Ghouzam                                   *
 *   Copyright (C) 2017   Guillaume Postic                                  *
 *   Copyright (C) 2017   Hubert Santuz                                     *
 *   guillaume.postic@univ-paris-diderot.fr                                 *
 *   jean-christophe.gelly@univ-paris-diderot.fr                            *
 *                                                                          *
 *   This program is free software: you can redistribute it and/or modify   *
 *   it under the terms of the GNU General Public License as published by   *
 *   the Free Software Foundation, either version 3 of the License, or      *
 *   (at your option) any later version.                                    *
 *                                                                          *
 *   This program is distributed in the hope that it will be useful,        *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *   GNU General Public License for more details.                           *
 *                                                                          *
 *   You should have received a copy of the GNU General Public License      *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                          *
 ***************************************************************************/

/* solvation.cpp */

#include "solvation.h"
#include <cstdlib>
#include <fstream> 
#include <sstream> 
#include <iostream> 

using namespace std;

#define NBR_DIST 34

Solvation::Solvation(const string &filename, float _w){
    weight=_w;
    ifstream sol_file(filename.c_str());
    string line="";
    istringstream iss;
    string mot,key;
    for(int i=0;i<20;i++){
        getline(sol_file,line);
        iss.str(line);
        iss>>key;
        for(int j=0;j<NBR_DIST;j++){
            iss >> mot;
            ener_residues[key].push_back(atof(mot.c_str()));
        }
        iss.clear();
    }
}


void Solvation::show(const string & key){
    cout<<key<<" :";
    for(int i=0; i<NBR_DIST;i++){
        cout<<" "<< this->ener_residues[key][i];
    }
    cout<<endl;
}


float Solvation::getEnergy(string key, int nbr){
    float energy=0;
    
    if(nbr>=NBR_DIST){
        energy=6.00;
    }else{
        energy=weight*this->ener_residues[key][nbr];
    }

    return energy;
}
