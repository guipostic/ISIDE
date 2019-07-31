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

/* parsetm.h */

#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>


/* SPLIT */
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}


/* SYSCALL */
// Returns the STDOUT as a string, including '\n'
std::string exec(char* cmd) {
    FILE* pipe = popen(cmd, "r");
    if (!pipe) return "ERROR";
    char buffer[128];
    std::string result = "";
    while(!feof(pipe)) {
        if(fgets(buffer, 128, pipe) != NULL)
                result += buffer;
    }
    pclose(pipe);
    return result;
}


/* Parse TMscore output */
std::vector<std::string> parseTM(char* command){
    std::vector<std::string> results;
    std::string syscallout = exec(command);
    std::vector<std::string> lines = split(syscallout, '\n');
    for(auto line : lines){
        if (line.substr(0,5) == "RMSD "
            || line.substr(0,8) == "TM-score"
            || line.substr(0,6) == "MaxSub"
            || line.substr(0,6) == "GDT-TS"
            || line.substr(0,6) == "GDT-HA"){
            std::vector<std::string> sides = split(line, '=');
            std::string value = sides[1];
            remove_if(value.begin(), value.end(), ::isspace);
            results.push_back(value);
        }
    }

    return results;
}
