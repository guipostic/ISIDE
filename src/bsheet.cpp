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

/* bsheet.cpp */

#include <cstdlib> 
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "bsheet.h"

using namespace std;

Bsheet::Bsheet(string name_para, string name_anti, float min_dist, float max_dist, float _nbr, float _w){
    weight=_w;
    dmin=min_dist;dmax=max_dist;nbr=_nbr;
    ifstream paraf(name_para.c_str());
    ifstream antif(name_anti.c_str());
    string line="";
    istringstream sline;
    float elem;
    
    for(int i=0;i<nbr;i++){
        getline(paraf,line);
        sline.str(line);
        ener_mat_para.push_back(vector<float>());

        for(int j=0;j<nbr;j++){
            sline >> elem;
            ener_mat_para[i].push_back(elem);
        }
        sline.clear();
    }

    for(int i=0;i<nbr;i++){
        getline(antif,line);
        sline.str(line);
        ener_mat_anti.push_back(vector<float>());

        for(int j=0;j<nbr;j++){
            sline >> elem;
            ener_mat_anti[i].push_back(elem);
        }
        sline.clear();
    }
}

float Bsheet::getEnergy(float dist, float angle){
    double ener=0;
    if(isnanf(angle)){
        cout<<"BUG_ANGLE"<<endl;
        return 0;
    }
    if(dist<dmin || dist>dmax){
        return 0;
    }
    if(angle<90){
        ener=weight*ener_mat_para[(int)(angle*(nbr/180.0))][(int)((dist-dmin)*(100.0/(dmax-dmin)))]; //distances allant de 2.5 à 7.0 ????????
    }else{
        ener=weight*ener_mat_anti[(int)(angle*(nbr/180.0))][(int)((dist-dmin)*(100.0/(dmax-dmin)))];
    }
    // DEBUG
    /*if(ener < -1000.0 || ener > 0.1){
        cout<< "Ener : "<<ener<<" D"<<(int)((dist-2.5)*(100.0/4.5))<<" : "<<dist<<" A"<<(int)(angle*(10.0/18.0))<<": "<<angle<<endl;
    }*/
    return ener;
}

