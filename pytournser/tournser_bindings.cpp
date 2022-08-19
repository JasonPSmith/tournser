#include <stdio.h>
#include <iostream>

#include "../tournser.cpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(pytourn, m) {

  m.doc() = "Python interface for tournser";

  m.def("run_tournser", [](std::vector<int>& vertices, std::vector<std::vector<value_t>>& edges,
                                char* filter, bool approx, uint32_t approx_val,
                                bool count_only, bool in_file, char* filename) {

    // Save std::cout status and disable
    auto cout_buff = std::cout.rdbuf();
    std::cout.rdbuf(nullptr);

    std::vector<size_t> cell_counts;
    py::dict output;
    delta_complex_t complex = delta_complex_t();

    //initialise approximate functionality
    size_t max_entries = std::numeric_limits<size_t>::max();
    if (approx) max_entries = approx_val;

    bool dist = false;
    if(in_file){
        complex = read_graph(filename,filter,dist);
    } else{
        complex = read_graph(vertices,edges,filter,dist);
    }
    complex.compute_oldest_cofaces();



    for(int i = 0; i < complex.cells.size(); i++){
        cell_counts.push_back(complex.cells[i].size());
    }
    
    output["cell_counts"] = cell_counts;

    //create tournser object and compute persistent homology
    if (!count_only){
        tournser thetournser = tournser(&complex,filename,max_entries,true);
        thetournser.compute_barcodes();
        output["finite_pairs"] = thetournser.finite_pairs;
        output["infinite_pairs"] = thetournser.infinite_pairs;
        output["bettis"] = thetournser.num_infinite_pairs();
    }
    
 // Re-enable again cout
    std::cout.rdbuf(cout_buff);

    return output;
  });
}
