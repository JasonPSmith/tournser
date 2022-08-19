#include <stdio.h>
#include <iostream>

#include "../tournser.cpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(pytournser, m) {

  m.doc() = "Python interface for tournser";

  m.def("run_tournser", [](std::vector<int>& vertices, std::vector<std::vector<value_t>>& edges,
                                char* filter, bool approx, uint32_t approx_val,
                                bool count_only, bool in_file, char* filename) {

    // Save std::cout status and disable
    auto cout_buff = std::cout.rdbuf();
    std::cout.rdbuf(nullptr);

    std::vector<size_t> cell_counts;
    delta_complex_t complex

    bool count_only = false;
    if (!count_only){
        std::ifstream f(outname);
        if (f.good()){
            std::cerr << "The output file already exists, aborting." << std::endl;
            exit(-1);
        } else { std::cout << "Printing output to : " << outname << std::endl; }
    }

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
    
    // Re-enable again cout
    std::cout.rdbuf(cout_buff);

    py::dict output;
    output["cell_counts"] = cell_counts;

    //create tournser object and compute persistent homology
    if (!count_only){
        std::ostringstream ss;
        tournser(&complex,outname,max_entries,true,ss).compute_barcodes();
        output["homology"] = ss.str();
    }
    

    return output;
  });
}
