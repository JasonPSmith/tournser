//Run with:
//./tournser in_file out_file
//plus any of the following optional inputs:
//--filtration     followed by one of: local, global, 3cycle, vertexMax, vertexSum, natural, natural-max
//--approximate    followed by an positive integer, skips any bars taking to long to compute, 100000 is a good value
//--print          followed by address of file to print into
//--print_dist     followed by true to print the distribution of filtrations
//--count_only     followed by true to skip homology computation, and no out_file is required
//-------------------------------------------------------------------------//
//BEGIN header

#define INDICATE_PROGRESS
#define PARALLEL_THREADS 8

#include <cassert>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <queue>
#include <sstream>
#include <unordered_map>
#include <set>
#include <functional>
#include <cmath>
#include <cstring>
#include <thread>

typedef float value_t;
typedef int64_t index_t;
typedef int16_t coefficient_t;
typedef index_t entry_t;
typedef std::pair<value_t, index_t> filtration_index_t;
typedef std::unordered_map<std::string, const char*> arg_t;

typedef std::deque<index_t> pivot_column_index_t;
const index_t INVALID_INDEX = std::numeric_limits<index_t>::max();

float string_to_float(std::string s) { return atof(s.c_str()); }

const index_t get_index(const entry_t& i) { return i; }
index_t get_coefficient(const entry_t& i) { return 1; }
entry_t make_entry(index_t _index, coefficient_t _value) { return entry_t(_index); }

const entry_t& get_entry(const entry_t& e) { return e; }

value_t get_filtration(const filtration_index_t& i) { return i.first; }
index_t get_index(const filtration_index_t& i) { return i.second; }

class filtration_entry_t : public std::pair<value_t, entry_t> {
public:
	filtration_entry_t() {}
	filtration_entry_t(const entry_t& e) : std::pair<value_t, entry_t>(0, e) {}
    filtration_entry_t(index_t _filtration, index_t _index)
        : std::pair<value_t, entry_t>(_filtration, entry_t(_index)) {}
	filtration_entry_t(value_t _filtration, index_t _index, coefficient_t _coefficient)
	    : std::pair<value_t, entry_t>(_filtration, make_entry(_index, _coefficient)) {}
	filtration_entry_t(const filtration_index_t& _filtration_index, coefficient_t _coefficient)
	    : std::pair<value_t, entry_t>(get_filtration(_filtration_index),
	                                 make_entry(get_index(_filtration_index), _coefficient)) {}
	filtration_entry_t(const filtration_index_t& _filtration_index)
	    : filtration_entry_t(_filtration_index, 1) {}
};

template <typename Heap> filtration_entry_t pop_pivot(Heap& column) {
    if (column.empty())
        return filtration_entry_t(-1);
    else {
        auto pivot = column.top();
        column.pop();
        while (!column.empty() &&
               get_index(column.top()) == get_index(pivot)) {
            column.pop();

            if (column.empty())
                return filtration_entry_t(-1);
            else {
                pivot = column.top();
                column.pop();
            }
        }
        return pivot;
    }
}

const entry_t& get_entry(const filtration_entry_t& p) { return p.second; }
const index_t get_index(const filtration_entry_t& p) { return get_index(get_entry(p)); }
const coefficient_t get_coefficient(const filtration_entry_t& p) { return get_coefficient(get_entry(p)); }
const value_t& get_filtration(const filtration_entry_t& p) { return p.first; }

template <typename Entry> struct greater_filtration_or_smaller_index {
	bool operator()(const Entry& a, const Entry& b) {
		return (get_filtration(a) > get_filtration(b)) ||
		       ((get_filtration(a) == get_filtration(b)) && (get_index(a) < get_index(b)));
	}
};

template <typename Entry> struct smaller_index {
	bool operator()(const Entry& a, const Entry& b) { return get_index(a) < get_index(b); }
};

std::string toBinary(int num, int bits) {
    std::vector<char> digits(bits);
    for (int i = 0; i < bits; ++i) {
        digits.push_back(num % 2 + '0');
        num >>= 1;
    }
    return std::string(digits.rbegin(), digits.rend());
}

void print_progress(int dim,int opt){
#ifdef INDICATE_PROGRESS
    if(opt == 0){ std::cout << "Constructing Dimension " << dim << std::endl; }
    else if(opt == 1){ std::cout << "Computing Dimension " << dim << std::endl; }
    else if(opt == 2){ std::cout << "Computing Homology :" << std::endl; }
    else if(opt == 3){ std::cout << "Tournaplex Fully Constructed : " << std::endl; }
    else if(opt == 4){ std::cout << "Constructing Tournaplex: " << std::endl; }
#endif
}

class filtered_union_find {
	std::vector<index_t> parent;
	std::vector<std::vector<index_t>> rank;
	const std::vector<value_t> filtration;

public:
	filtered_union_find(const std::vector<value_t>& _filtration)
	    : rank(_filtration.size()), filtration(_filtration), parent(_filtration.size()) {
		for (index_t i = 0; i < _filtration.size(); ++i){
			parent[i] = i;
			rank[i] = std::vector<index_t>{i};
		}
	}
	index_t find(index_t x) {
		return parent[x];
	}
	value_t link(index_t x, index_t y) {
		x = find(x);
		y = find(y);
		if (x == y) return -1;
		if (filtration[x] < filtration[y] || (filtration[x] == filtration[y] && rank[x].size() > rank[y].size())){
			for(auto i : rank[y]) parent[i] = x;
			rank[x].insert(rank[x].end(),rank[y].begin(),rank[y].end());
			return filtration[y];
		} else {
			for(auto i : rank[x]) parent[i] = y;
			rank[y].insert(rank[y].end(),rank[x].begin(),rank[x].end());
			return filtration[x];
		}
	}
	value_t filter(index_t x){
		return filtration[find(x)];
	}
};

template <typename ValueType> class compressed_sparse_matrix {
	std::deque<size_t> bounds;
	std::deque<ValueType> entries;

public:
	size_t size() const { return bounds.size(); }

	void clear() {
		bounds.clear();
		bounds.shrink_to_fit();
		entries.clear();
		entries.shrink_to_fit();
	}

	typename std::deque<ValueType>::const_iterator cbegin(size_t index) const {
		assert(index < size());
		return index == 0 ? entries.cbegin() : entries.cbegin() + bounds[index - 1];
	}

	typename std::deque<ValueType>::const_iterator cend(size_t index) const {
		assert(index < size());
		return entries.cbegin() + bounds[index];
	}

	template <typename Iterator> void append_column(Iterator begin, Iterator end) {
		for (Iterator it = begin; it != end; ++it) { entries.push_back(*it); }
		bounds.push_back(entries.size());
	}

	void append_column() { bounds.push_back(entries.size()); }

	void push_back(ValueType e) {
		assert(0 < size());
		entries.push_back(e);
		++bounds.back();
	}

	void pop_back() {
		assert(0 < size());
		entries.pop_back();
		--bounds.back();
	}

	template <typename Collection> void append_column(const Collection collection) {
		append_column(collection.cbegin(), collection.cend());
	}
};
int binom(int n, int k) {
   if (k == 0 || k == n){ return 1; }
   return binom(n - 1, k - 1) + binom(n - 1, k);
}

//END header
//-------------------------------------------------------------------------//
//BEGIN delta_complex

//TODO Make most of the elements private for all three classes


class delta_complex_t;

//removes trailing white space
static inline void trim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

template <typename t>
std::vector<t> split(std::string& s, char delim, const std::function<t(std::string)>& transform) {
	trim(s);
	std::vector<t> elems;
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)){ elems.push_back(transform(item)); }
	return elems;
}

class delta_complex_cell_t {
public:
	//Stored as the vertices and an orientation of the edges,
	//With orientations stored as a lower triangular matrix (via a vector of vectors) in the following way
	//orientation[a-1][b] == 0 iff edge orientated vertices[a]->vertices[b] (always assuming a>b)
	//We also store the boundary and coboundary as their indices in cells[dim-1] and cells[dim+1], resp.
    std::vector<index_t> coboundary;
	std::vector<index_t> boundary;
	std::vector<index_t> vertices;
	std::vector<std::vector<index_t>> orientation;
	index_t location;
    index_t cli_loc;
	value_t filt;
	index_t oldest_coface;
	//If cell is a vertex then we also store the neighbours in the original graph, this is used when constructing the tournaplex
	//neighbours is a set of all neighbours and in_or_out indicates whether they are a in-neighbour or out-neighbour
	//for each neighbour u returns in_or_out[u] returns 0 if current->u is an edge, 1 if u->current is edge and 2 if both are edges
	std::set<index_t> neighbours;
	std::unordered_map<index_t,int> in_or_out;

    //initialise class
	delta_complex_cell_t(index_t _vertex, value_t _filt) : filt(_filt),
						oldest_coface(-1), vertices(std::vector<index_t>{_vertex}), location(_vertex) {}
    delta_complex_cell_t(std::vector<index_t> _vertices, std::vector<std::vector<index_t>> _orientation, value_t _filt, index_t _loc) : filt(_filt),
                        oldest_coface(-1), orientation(_orientation), vertices(_vertices), location(_loc), cli_loc(-1) {}
	delta_complex_cell_t(std::vector<index_t> _vertices, std::vector<std::vector<index_t>> _orientation, value_t _filt, index_t _loc, index_t _cli_loc) : filt(_filt),
						oldest_coface(-1), orientation(_orientation), vertices(_vertices), location(_loc), cli_loc(_cli_loc) {}

	bool same_orientation(std::vector<std::vector<index_t>> in){
       return orientation == in;
    }
	index_t dimension(){
		return vertices.size()-1;
	}
	void add_neighbour(index_t to_add, int b){
		neighbours.insert(to_add);
		if(in_or_out.size() < neighbours.size()) in_or_out[to_add] = b;
		else in_or_out[to_add] = 2;
	}
	value_t filtration(){
		return filt;
	}
	void print_cell(std::ofstream& outstream){
		outstream << "Vrts = ";
		for(auto v : vertices) outstream << v << " ";
		outstream << "| bdry = ";
		for(auto b : boundary) outstream << b << " ";
		outstream << "| cbdry = ";
		for(auto c : coboundary) outstream << c << " ";
		outstream << "| ort = ";
		for(auto d : orientation){
			for(auto e : d){
				outstream << e << " ";
			}
			outstream << ": ";
		}
		outstream << "| filt = " << filt << " | loc = " << location << std::endl;
	}
	void compute_boundary(delta_complex_t* complex, size_t thread, std::vector<std::vector<std::pair<index_t, index_t>>>& new_cell_cofaces);
	void compute_oldest_coface(delta_complex_t* complex);
	void assign_filtration(const char* filter,delta_complex_t* complex);
}; //END delta_complex_cell_t

class clique_t {
//This class is used when constructing the tournaplex but then cleared
//It stores the vertices that make a clique and the index of all faces with that vertex set
//via their indices in cells[dim]
//We also store the map cofaces which given any vertex u large than max(vertices)
//cofaces[u] returns the index of the clique {vertices} union {u} in cliques[dim+1]
public:
	std::vector<index_t> faces;
	std::vector<index_t> vertices;
	index_t location;

	clique_t(){}
	clique_t(index_t _vertex) : location(_vertex) {vertices.push_back(_vertex);faces.push_back(_vertex);}
	clique_t(std::vector<index_t> _vertices, index_t _loc)  : vertices(_vertices), location(_loc) {}

	index_t dimension(){
		return vertices.size()-1;
	}
	bool same_vertices(std::vector<index_t> in){
		return in == vertices;
	}
	void add_face(index_t cell){
		faces.push_back(cell);
	}
	void add_coface(index_t a,index_t b){
		cofaces.insert(std::make_pair(a,b));
	}
	index_t get_coface(index_t in){
		auto loc = cofaces.find(in);
		if(loc == cofaces.end()) return -1;
		return loc->second;
	}
	void print_cofaces(){
		for(auto d : cofaces){ std::cout << d.first<<","<<d.second << " "; }
	}
private:
	std::unordered_map<index_t,index_t> cofaces;
};

class delta_complex_t {
public:
    //stores all cells as a vector of vectors, each vector being all cells of that dimension
	//Also store all cliques in a similar fashion, but this is cleared after the tournaplex is constructed
	std::vector<std::vector<delta_complex_cell_t>> cells;
    std::vector<std::vector<clique_t>> cliques;


	//creates a delta complex with with |vertices| number of vertices
	//and each vertex i has filtration value given by vertices[i]
	//also create a clique of size one for each vertex
	delta_complex_t(){};
	delta_complex_t(std::vector<int> vertices){
		add_dimension();
		for( int i = 0; i < vertices.size(); i++ ){
            cells[0].push_back(delta_complex_cell_t(i,vertices[i]));
			cliques[0].push_back(clique_t(i));
        }
	};

	void print_complex(std::string filename){
        std::ofstream outstream;
        outstream.open(filename);
		for(int i = 0; i < cells.size(); i++){
			outstream << "Dim " << i << std::endl;
			for(int j = 0; j < cells[i].size(); j++){
				cells[i][j].print_cell(outstream);
			}
		}
        outstream.close();
	}

	void add_edge(index_t i,index_t j,value_t val){
		//NOTE: We assume edges are added in increasing lexicographic order,
		//this avoids always checking if relevant 2-clique already in list
		//This could be changed and it would make the code easier to read,
		//But I think it is slightly faster this way, though probably a negligible amount
		index_t loc = cells[1].size();
		std::vector<std::vector<index_t>> new_orient;
		new_orient.push_back(std::vector<index_t>());
		if( i < j ){
			//When i < j create a new clique and a new face, with edge orientation 1
			//add the new face to the new clique and also add the cofaces to the clique
			new_orient[0].push_back(1);
			cells[1].push_back(delta_complex_cell_t(std::vector<index_t>{i,j},new_orient,val,loc));
			index_t cli_loc = cliques[1].size();
			cliques[1].push_back(clique_t(std::vector<index_t>{i,j},cli_loc));
			cliques[1].back().faces.push_back(loc);
			cliques[0][i].add_coface(j,cli_loc);
		} else{
			//When i > j  create a new face with orientation 0
			//then check if the clique already exists, if it does not then create it
			//then add the new face the the clique
			new_orient[0].push_back(0);
			cells[1].push_back(delta_complex_cell_t(std::vector<index_t>{j,i},new_orient,val,loc));
			index_t the_cli = get_clique(std::vector<index_t>{j,i});
			if(the_cli == -1){
				index_t cli_loc = cliques[1].size();
				cliques[1].push_back(clique_t(std::vector<index_t>{j,i},cli_loc));
				cliques[0][j].add_coface(i,cli_loc);
				the_cli = cli_loc;
			}
			cliques[1][the_cli].faces.push_back(loc);
		}
		//Add the relevent neighbour and coboundary information to the vertices
		cells[0][i].add_neighbour(j,0);
		cells[0][j].add_neighbour(i,1);
		cells[0][i].coboundary.push_back(loc);
		cells[0][j].coboundary.push_back(loc);
		//Add the vertices of faces of the new created face
		cells[1].back().boundary.push_back(i);
		cells[1].back().boundary.push_back(j);
	}
	void add_dimension(){
		//increase dimension of complex by adding new layer to cells and cliques
		cliques.push_back( std::vector<clique_t>() );
		cells.push_back( std::vector<delta_complex_cell_t>() );
	}
	index_t get_clique(std::vector<index_t> nodes){
		//Returns the index in cliques[|nodes|] of the clique with vertex set nodes, or returns -1 if no such clique exists.
		//Let nodes be v_1,v_2,.., note they are assumed to be inputted in increasing order
		//It does this by looking at vertex face v_0 and asking that where the face v_0,v_1 is
		//and then asking that where the face v_0,v_1,v_2 is etc until we find the face
		clique_t* cli = &cliques[0][nodes[0]];
		for(int i = 1; i < nodes.size(); i++){
			index_t x = cli->get_coface(nodes[i]);
			if(x == -1) return -1;
			cli = &cliques[i][x];
		}
		return cli->location;
	}
	void build_tournaplex(const char* filter){
		//This constructs the tournaplex from its 1-skeleton
		int dim = 2;
		bool go = true;
		while(go){
            print_progress(dim,0);
			add_dimension();

            std::thread t[PARALLEL_THREADS - 1];
            const size_t cliques_per_thread = cliques[dim-1].size() / PARALLEL_THREADS;
            std::vector<std::vector<delta_complex_cell_t>> new_cells(PARALLEL_THREADS,std::vector<delta_complex_cell_t>());
            std::vector<std::vector<clique_t>> new_cliques(PARALLEL_THREADS,std::vector<clique_t>());
            std::vector<std::vector<std::tuple<index_t, index_t, index_t>>> new_clique_cofaces(PARALLEL_THREADS,std::vector<std::tuple<index_t, index_t, index_t>>());
            std::vector<std::vector<std::pair<index_t, index_t>>> new_cell_cofaces(PARALLEL_THREADS,std::vector<std::pair<index_t, index_t>>());

			//For each clique of dimension one less, pick a vertex loop through its neighbours
			//for each neighbour we check if it is  neighbour of all vertices in the clique, if it is then add the relevant faces
            for (size_t index = 0; index < PARALLEL_THREADS - 1; ++index){
    			t[index] = std::thread(&delta_complex_t::if_common_create_all, this, index*cliques_per_thread, (index+1)*cliques_per_thread,
                                        dim, filter, index, std::ref(new_cells), std::ref(new_cliques), std::ref(new_clique_cofaces), std::ref(new_cell_cofaces));
            }
            if_common_create_all((PARALLEL_THREADS-1)*cliques_per_thread,cliques[dim-1].size(), dim, filter, PARALLEL_THREADS-1,
                                                                            new_cells, new_cliques, new_clique_cofaces, new_cell_cofaces);

            //wait until threads stop
            for (size_t i = 0; i < PARALLEL_THREADS - 1; ++i) t[i].join();

            //We combine the outputs from each threads
            int cell_start = 0;
            int clique_start = 0;
            for (size_t index = 0; index < PARALLEL_THREADS; ++index){
                for(int i = 0; i < new_cliques[index].size(); i++){
                    cliques.back().push_back(new_cliques[index][i]);
                    cliques.back().back().location = new_cliques[index][i].location + clique_start;
                    for(auto f : cliques.back().back().faces){
                        f += cell_start;
                    }
                    cliques[dim-1][std::get<0>(new_clique_cofaces[index][i])].add_coface(std::get<1>(new_clique_cofaces[index][i]),std::get<2>(new_clique_cofaces[index][i])+clique_start);
                }
                for(auto c : new_cells[index]){
                    cells.back().push_back(c);
                    cells.back().back().location = c.location+cell_start;
                    cliques.back()[c.cli_loc+clique_start].faces.push_back(c.location+cell_start);
                }
                for(auto c : new_cell_cofaces[index]){
                    cells[dim-1][c.first].coboundary.push_back(c.second+cell_start);
                }

                cell_start = cells.back().size();
                clique_start = cliques.back().size();
            }

			//Once we reach a pass where no new faces were create we delete the empty dimension
			//clear cliques from memory and finish
			if(cells[dim].size() == 0){
				go = false;
				cells.pop_back();
				cliques.clear();
			} else {dim++;}
		}

	}
    void if_common_create_all(size_t start, size_t end, int dim, const char* filter, size_t thread,
         std::vector<std::vector<delta_complex_cell_t>>& new_cells, std::vector<std::vector<clique_t>>& new_cliques,
        std::vector<std::vector<std::tuple<index_t, index_t, index_t>>>& new_clique_cofaces,std::vector<std::vector<std::pair<index_t, index_t>>>& new_cell_cofaces){
        for(int c = start; c < end; c++){
            for(auto v : cells[0][cliques[dim-1][c].vertices[0]].neighbours){
                if_common_create(c,v,dim,filter,thread,new_cells,new_cliques,new_clique_cofaces,new_cell_cofaces);
            }
        }
    }
	void if_common_create(index_t cl, index_t v, index_t dim, const char* filter, size_t thread,
        std::vector<std::vector<delta_complex_cell_t>>& new_cells, std::vector<std::vector<clique_t>>& new_cliques,
        std::vector<std::vector<std::tuple<index_t, index_t, index_t>>>& new_clique_cofaces,std::vector<std::vector<std::pair<index_t, index_t>>>& new_cell_cofaces ){
		//This checks if the vertex v is a neighbour of all vertices in the clique cliques[dim-1][cl]
		//If it is we create a new clique and then compute all the orientations of the edges and create face for each

		//If v is smaller than the largest vertex in the clique we can ignore it,
		//because such a clique would have been constructed earlier
		//as we construct them by adding larger valued vertices at each step
		if(v <= cliques[dim-1][cl].vertices.back()) return;

		//For each vertex u in the clique check if it is a neighbour of v,
		//if it is the record the orientation of the edge(s) between u and v
		std::vector<int> orient;
		orient.push_back(cells[0][v].in_or_out.find(cliques[dim-1][cl].vertices[0])->second);
		for(int i = 1; i < cliques[dim-1][cl].vertices.size(); i++){
			auto it = cells[0][v].in_or_out.find(cliques[dim-1][cl].vertices[i]);
			if(it == cells[0][v].in_or_out.end()) return;
			orient.push_back(it->second);
		}

		//If we reach this point we know that v is a neighbour of all vertices in the clique
		//so we create a new clique by adding v to the previous vertex set
		//Store the cofaces we need to add to old clique, we do this later as need to know index once parrallel threads have combined
		std::vector<index_t> new_vertices;
		for(auto j : cliques[dim-1][cl].vertices) new_vertices.push_back(j);
		new_vertices.push_back(v);
		index_t cli_loc = new_cliques[thread].size();
		new_cliques[thread].push_back(clique_t(new_vertices,cli_loc));
        new_clique_cofaces[thread].push_back(std::make_tuple(cl,v,cli_loc));

		//compute number of bidirectional edges in new clique
		int num_bidirect = 0;
		for(int k = 0; k < orient.size(); k++){
			if(orient[k] == 2) num_bidirect++;
		}

		//take each number between 0 and the number of bidirection edges and represent it in binary_function
		//we use this to create all orientations of the edges
		for(int i = 0; i < pow(2,num_bidirect); i++){
			//create a new orientation matrix by taking the orientation matrix of the old clique_lookup
			//and add an new bottom row, in this new bottom row we orientate any unidirectional edges
			//in there unique way and then use the binary number to decide how the bidirectional  edges are orientated
			std::vector<index_t> top_row;
			std::string bin = toBinary(i,num_bidirect);
			int count = 0;
			for(int j = 0 ; j < orient.size(); j++){
				if(orient[j] == 2){
					if(bin[count] == '0'){ top_row.push_back(0);
					} else{ top_row.push_back(1); }
					count++;
				} else{
					top_row.push_back(orient[j]);
				}
			}
			//Here we add our new bottom row to each of the orientations of original clique
			for(auto f : cliques[dim-1][cl].faces){
				std::vector<std::vector<index_t>> new_orients;
				for(auto t : cells[dim-1][f].orientation) new_orients.push_back(t);
				new_orients.push_back(top_row);
				//We then create a new face with this orientation matrix
				//add out new face the new clique, and compute the boundary of the new face
				index_t loc = new_cells[thread].size();
				new_cells[thread].push_back(delta_complex_cell_t(new_vertices,new_orients,0,loc,cli_loc));
				new_cells[thread].back().compute_boundary(this,thread,new_cell_cofaces);
				new_cells[thread].back().assign_filtration(filter,this);
			}
		}
    }
    index_t number_of_cells(index_t dimension) const {
        if ( dimension >= cells.size() ) { return 0; }
        return cells[dimension].size();
    }
    bool is_top_dimension(index_t dimension){
        return cells.size() == dimension-1;
    }
    index_t top_dimension(){
        return cells.size()-1;
    }
    value_t filtration(index_t dimension,index_t index){
        return cells[dimension][index].filtration();
    }
	std::vector<value_t> vertex_filtration(){
		std::vector<value_t> out;
		for( auto c : cells[0] ){
			out.push_back(c.filtration());
		}
        return out;
    }
    delta_complex_cell_t* get(index_t dimension,index_t index){
        return &cells[dimension][index];
    }
	void compute_oldest_cofaces(){
		for (auto p : cells){
			for(auto q : p){
				q.compute_oldest_coface(this);
			}
		}
	}
	void initialise_filtration(const char* filter){
		for(int i = 0; i < cells[0].size(); i++){
			cells[0][i].assign_filtration(filter,this);
		}
		for(int i = 0; i < cells[1].size(); i++){
			cells[1][i].assign_filtration(filter,this);
		}
	}
//END delta_complex_t
};

void delta_complex_cell_t::compute_oldest_coface(delta_complex_t* complex){
	value_t oldest = -1;
	for( auto c : coboundary ){
		value_t f = complex->get(dimension()+1,c)->filtration();
		if( f > oldest ){
			oldest = f;
			oldest_coface = c;
		}
	}
}
void delta_complex_cell_t::compute_boundary(delta_complex_t* complex, size_t thread, std::vector<std::vector<std::pair<index_t, index_t>>>& new_cell_cofaces){
	//For each vertex delete it from the face, by removing it from the vertex set and
	//the relevant rows and columns of the orientation matrix
	//find this face in the complex and add its index to the boundary of the cell,
	//also add this cell to the coboundary of the found cell
	for(int i = 0; i < vertices.size(); i++){
		std::vector<index_t> temp;
		std::vector<std::vector<index_t>> tempOrient;
		for(int j = 0; j < vertices.size(); j++){
			if(j != i){
				temp.push_back(vertices[j]);
				if(j < i && j > 0){ tempOrient.push_back(orientation[j-1]);
				} else if(j > i && j > 1){
					tempOrient.push_back(std::vector<index_t>());
					for(int k = 0; k < orientation[j-1].size(); k++){
						if(k != i){
							tempOrient.back().push_back(orientation[j-1][k]);
						}
					}
				}
			}
		}
        bool found = false;
        int f = 0;
        while(!found && f < complex->cliques[dimension()-1][complex->get_clique(temp)].faces.size()){
			if(complex->cells[dimension()-1][complex->cliques[dimension()-1][complex->get_clique(temp)].faces[f]].same_orientation(tempOrient)){
				boundary.push_back(complex->cells[dimension()-1][complex->cliques[dimension()-1][complex->get_clique(temp)].faces[f]].location);
                new_cell_cofaces[thread].push_back(std::make_pair(complex->cliques[dimension()-1][complex->get_clique(temp)].faces[f],location));
				found = true;
			}
            f++;
		}
	}
}
void delta_complex_cell_t::assign_filtration(const char* filter,delta_complex_t* complex){
	if(strcmp(filter,"local") == 0){
		if(dimension() == 0){ filt = 0; }
		else if(dimension() == 1){ filt = 2; }
		else{
			std::vector<index_t> degrees(vertices.size());
			for(int i = 0; i < orientation.size(); i++){
				for(int j = 0; j < orientation[i].size(); j++){
				    if(orientation[i][j] == 0){
						degrees[i+1]++;
						degrees[j]--;
					}
					else{
					   degrees[i+1]--;
					   degrees[j]++;
					}
				}
			}
			int x = 0;
			for(auto i : degrees){ x += pow(i,2); }
			int n = dimension();
			x += ((n+1)*n*(n-1)/3);
			filt = x;
		}
	} else if(strcmp(filter,"global") == 0){
		if(dimension() == 0){
			value_t m = 0;
			for( auto i : in_or_out ){
				if( i.second == 0 ) { m += 1; }
				else if( i.second == 1 ) { m -= 1; }
			}
			filt = pow(m,2);
		} else if(dimension() > 0){
			value_t m = complex->cells[0][vertices[0]].filtration();
			for(int i = 1; i < vertices.size(); i++){
				m += complex->cells[0][vertices[i]].filtration();
			}
			filt = m;
		}
	}else if(strcmp(filter,"vertexMax") == 0){
		if(dimension() > 0){
			value_t m = complex->cells[0][vertices[0]].filtration();
			for(int i = 1; i < vertices.size(); i++){
				m = std::max(m,complex->cells[0][vertices[i]].filtration());
			}
			filt = m;
		}
	}else if(strcmp(filter,"edgeMax") == 0){
		if(dimension() > 1){
			value_t m = complex->cells[1][vertices[0]].filtration();
			for(int i = 1; i < vertices.size(); i++){
				m = std::max(m,complex->cells[1][vertices[i]].filtration());
			}
			filt = m;
		}
	}else if(strcmp(filter,"vertexSum") == 0){
		if(dimension() > 0){
			value_t m = complex->cells[0][vertices[0]].filtration();
			for(int i = 1; i < vertices.size(); i++){
				m += complex->cells[0][vertices[i]].filtration();
			}
			filt = m;
		}
	}else if(strcmp(filter,"3cycle") == 0){
		if(dimension() == 0){ filt = 0; }
		else if(dimension() == 1){ filt = 0; }
		else{
            this->assign_filtration("local",complex);
            int n = dimension()+1;
			filt = (binom(n,3)*((4*n-2)/(n-2)) - filt)/8;
		}
    }else if(strcmp(filter,"natural") == 0){
        if(dimension() == 0){ filt = 0; }
		else if(dimension() == 1){ filt = 2; }
		else{
			std::vector<index_t> degrees(vertices.size());
			for(int i = 0; i < orientation.size(); i++){
				for(int j = 0; j < orientation[i].size(); j++){
				    if(orientation[i][j] == 0){
						degrees[i+1]++;
						degrees[j]--;
					}
					else{
					   degrees[i+1]--;
					   degrees[j]++;
					}
				}
			}
			int x = 0;
			for(auto i : degrees){ x += pow(i,2); }
			value_t max_face = 0;
            for(int i = 1; i < boundary.size(); i++){
                max_face = std::max(max_face,complex->cells[dimension()-1][boundary[i]].filtration());
            }
			x += max_face;
			filt = x;
		}
    }else if(strcmp(filter,"natural-max") == 0){
        if(dimension() == 0){ filt = 0; }
        else if(dimension() == 1){ filt = 2; }
        else{
            std::vector<index_t> degrees(vertices.size());
            for(int i = 0; i < orientation.size(); i++){
                for(int j = 0; j < orientation[i].size(); j++){
                    if(orientation[i][j] == 0){
                        degrees[i+1]++;
                        degrees[j]--;
                    }
                    else{
                       degrees[i+1]--;
                       degrees[j]++;
                    }
                }
            }
            int x = 0;
            for(auto i : degrees){ x += pow(i,2); }
            value_t max_face = 0;
            for(int i = 1; i < boundary.size(); i++){
                max_face = std::max(max_face,complex->cells[dimension()-1][boundary[i]].filtration());
            }
            filt = std::max(x,(int)max_face);
        }
	} else{
		std::cout << "Error no matching filtration found, aborting." << std::endl;
		exit(-1);
	}
}

class simplex_coboundary_enumerator {
private:
	index_t  idx_above, dim;
	delta_complex_cell_t* simplex;
	delta_complex_t* complex;

public:
	simplex_coboundary_enumerator(const filtration_entry_t _simplex, index_t _dim,
	                              delta_complex_t* _complex)
	    : idx_above(0), simplex(_complex->get(_dim,_simplex.second)), dim(_dim), complex(_complex) {}

	bool has_next() {
		return idx_above < simplex->coboundary.size();
	}

	filtration_entry_t next() {
        idx_above++;
        return filtration_entry_t(complex->get(dim+1,simplex->coboundary[idx_above-1])->filtration(), simplex->coboundary[idx_above-1]);
	}
};

struct greater_filtration_or_better_pivot_or_smaller_index {
	greater_filtration_or_better_pivot_or_smaller_index(delta_complex_t* _complex, index_t _dimension) : complex(_complex), dimension(_dimension) {}
	bool operator()(filtration_index_t a, filtration_index_t b) const {
		// First order by the filtration value
		if (get_filtration(a) > get_filtration(b)) return true;
		if (get_filtration(a) < get_filtration(b)) return false;

		auto ta = get_coboundary_size_and_gap_and_pivot(get_index(a));
		auto tb = get_coboundary_size_and_gap_and_pivot(get_index(b));

		// Then the number of non-trivial coboundary entries
		if (std::get<0>(ta) < std::get<0>(tb)) return true;
		if (std::get<0>(ta) > std::get<0>(tb)) return false;

		// Then order by the better pivoting
		if (std::get<2>(ta) < std::get<2>(tb)) return true;
		if (std::get<2>(ta) > std::get<2>(tb)) return false;

		if (std::get<1>(ta) > std::get<1>(tb)) return true;
		if (std::get<1>(ta) < std::get<1>(tb)) return false;

		// Finally, order by their indices
		return get_index(a) < get_index(b);
	}

private:
	delta_complex_t* complex;
	index_t dimension;

	// A column is considered to be a better pivot if the jump from pivot to the next
	// non-trivial element is as big as possible. This prevents accidentally inserting
	// non-trivial elements just below the pivot, which sometimes creates very long
	// reduction chains.
	// The second sort criterium is for it to be small because the small pivots will be
	// used the most.
	std::tuple<size_t, size_t, index_t> get_coboundary_size_and_gap_and_pivot(index_t a) const {
		// Look at the first two gaps of the pivot and the next element
		index_t pivot = 0;
		size_t gap_after_pivot = 0;
		simplex_coboundary_enumerator iterator(a,dimension,complex);
		size_t coboundary_size = 0;
		while (iterator.has_next()) {
			coboundary_size++;
			index_t next_index = get_index(iterator.next().second);
			if (next_index > pivot) {
				gap_after_pivot = next_index - pivot;
				pivot = next_index;
			}
		}

		return std::make_tuple(coboundary_size, gap_after_pivot, pivot);
	}
};

// This class is just an ordinary priority queue, but once the
// queue gets too long (because a lot of faces are inserted multiple
// times) it starts collecting the coefficients and only inserting each
// new face once
template <class Container, class Comparator>
class priority_queue_t : public std::priority_queue<filtration_entry_t, Container, Comparator> {
	std::unordered_map<index_t, bool> coefficients;
	static const filtration_entry_t dummy;
	bool use_dense_version = false;
	size_t dense_threshold;

public:
	priority_queue_t(size_t _dense_threshold)
	    : dense_threshold(_dense_threshold) {}

	void push(const filtration_entry_t& value) {
		if (use_dense_version) {
			// If we already have this value: update the count and don't push it again
			auto p = coefficients.find(get_index(value));
			if (p != coefficients.end()) {
				p->second = !p->second;
				return;
			}
		}

		std::priority_queue<filtration_entry_t, Container, Comparator>::push(value);

		if (use_dense_version) coefficients.insert(std::make_pair(get_index(value), get_coefficient(value)));

		if (!use_dense_version &&
		    std::priority_queue<filtration_entry_t, Container, Comparator>::size() >= dense_threshold)
			use_dense_version = true;
	}

	void pop() {
		// Don't use this, only allow get_pivot
		throw std::exception();
	}

	filtration_entry_t pop_pivot() {
		remove_trivial_coefficient_entries();
		if (std::priority_queue<filtration_entry_t, Container, Comparator>::empty())
			return filtration_entry_t(-1);
		else {
			auto pivot = get_top();
			safe_pop();
			while (!std::priority_queue<filtration_entry_t, Container, Comparator>::empty() &&
			       get_index(std::priority_queue<filtration_entry_t, Container, Comparator>::top()) ==
			           get_index(pivot)) {
				safe_pop();
				remove_trivial_coefficient_entries();

				if (std::priority_queue<filtration_entry_t, Container, Comparator>::empty())
					return filtration_entry_t(-1);
				else {
					pivot = get_top();
					safe_pop();
				}
			}
			return pivot;
		}
	}

	filtration_entry_t get_pivot() {
		filtration_entry_t result = pop_pivot();
		if (get_index(result) != -1) { push(result); }
		return result;
	}

private:
	inline filtration_entry_t get_top() {
		auto pivot = std::priority_queue<filtration_entry_t, Container, Comparator>::top();
		return pivot;
	}

	inline void safe_pop() {
		if (use_dense_version) {
			auto e =
			    coefficients.find(get_index(std::priority_queue<filtration_entry_t, Container, Comparator>::top()));
			if (e != coefficients.end()) coefficients.erase(e);
		}
		std::priority_queue<filtration_entry_t, Container, Comparator>::pop();
	}

	inline void remove_trivial_coefficient_entries() {
		if (use_dense_version) {
			auto p = coefficients.find(get_index(std::priority_queue<filtration_entry_t, Container, Comparator>::top()));
			while (p != coefficients.end() && p->second == false) {
				coefficients.erase(p);
				std::priority_queue<filtration_entry_t, Container, Comparator>::pop();
				p = coefficients.find(get_index(std::priority_queue<filtration_entry_t, Container, Comparator>::top()));
			}
		}
	}
};
template <class Container, class Comparator>
const filtration_entry_t priority_queue_t<Container, Comparator>::dummy(filtration_entry_t(0.0, -1));



//END delta_complex
//-------------------------------------------------------------------------//
//BEGIN read_graph

void build_count_and_filter(delta_complex_t& complex, const char* filter, bool dist){
	complex.initialise_filtration(filter);
	complex.build_tournaplex(filter);

    print_progress(0,3);

	std::cout << "Simplex Count = ";
    for(int i = 0; i < complex.cells.size(); i++){
        std::cout << complex.cells[i].size() << " ";
    }
    std::cout << std::endl;

    //Compute and print the distribution of directionality
    if(dist){
        std::vector<std::unordered_map<entry_t,size_t>> direct_dist;
        std::cout << "Filtration Distribution (value,count): " << std::endl;
        for(int i = 0; i < complex.cells.size(); i++){
            direct_dist.push_back( std::unordered_map<entry_t,size_t>() );
            for(int j = 0; j < complex.cells[i].size(); j++){
                auto it = direct_dist[i].find(complex.cells[i][j].filtration());
                if(it != direct_dist[i].end()){ it->second = it->second+1; }
                else{ direct_dist[i].insert ( std::pair<entry_t,size_t>(complex.cells[i][j].filtration(),1) ); }
            }
            std::cout << "dim " << i << " : ";
            for(std::unordered_map<entry_t,size_t>::iterator it = direct_dist[i].begin();
                                                            it != direct_dist[i].end(); ++it) {
                std::cout << "("<< it->first <<","<< it->second <<") ";
            }
            std::cout << std::endl;
        }
    }
    std::cout << std::endl;
}

delta_complex_t read_graph(const std::string filename, const char* filter, bool dist){
	std::ifstream infile;
    infile.open(filename);
	std::string line;

	std::getline(infile,line);
	std::getline(infile,line);
    print_progress(0,0);

	std::vector<int> vertices = split<int>(line, ' ', string_to_float);
	delta_complex_t complex = delta_complex_t(vertices);
	complex.add_dimension();
	std::getline(infile,line);

    print_progress(1,0);
    int prev1 = std::numeric_limits<int>::min();
    int prev2 = std::numeric_limits<int>::min();
    while (std::getline(infile,line)){
        if (line.length() == 0) continue;
        else {
        	std::vector<int> faces = split<int>(line, ' ', string_to_float);
			value_t val = 0;
			if(faces.size() == 3) val = faces[2];
            if(faces[0] <= prev1 && faces[1] <= prev2){
                std::cout << "Error: edges are not in lexicographic order or added same edge twice." << std::endl;
                exit(-1);
            }
			complex.add_edge(faces[0],faces[1],val);
            prev1 = faces[0];
            prev2 = faces[1];
        }
    }
	build_count_and_filter(complex, filter, dist);

    return complex;

}

delta_complex_t read_graph(std::vector<int>& vertices, std::vector<std::vector<value_t>>& edges, const char* filter, bool dist){
	delta_complex_t complex = delta_complex_t(vertices);
	complex.add_dimension();

    print_progress(1,0);
    int prev1 = std::numeric_limits<int>::min();
    int prev2 = std::numeric_limits<int>::min();
    for (auto i : edges){
 		value_t val = 0;
		if(i.size() == 3) val = i[2];
        if(i[0] <= prev1 && i[1] <= prev2){
            std::cout << "Error: edges are not in lexicographic order or added same edge twice." << std::endl;
            exit(-1);
        }
		complex.add_edge(i[0],i[1],val);
        prev1 = i[0];
        prev2 = i[1];
    }
    build_count_and_filter(complex, filter, dist);

    return complex;
}

//END read_graph
//-------------------------------------------------------------------------//
//BEGIN tournser




class tournser {
	delta_complex_t* complex;
	std::ofstream outfile;
	std::ostringstream* ss;
	index_t n, dim_max;
	mutable std::vector<filtration_entry_t> coface_entries;
	size_t max_entries;
	std::vector<size_t> skipped_entries;
	bool python;

public:
    std::vector<std::vector<std::pair<value_t,value_t>>> finite_pairs;
	std::vector<std::vector<value_t>> infinite_pairs;

	tournser(	delta_complex_t* _complex, char* _outname, size_t _max_entries, bool _python)
	    : complex(_complex), n(complex->number_of_cells(0)),
	      dim_max(complex->top_dimension()),
		  max_entries(_max_entries),
		  python(_python){
				if(!python) {outfile = std::ofstream(_outname);}
				skipped_entries.assign(dim_max+1,0);
				infinite_pairs.resize(dim_max+1);
				finite_pairs.resize(dim_max+1);
			  }

	value_t compute_filtration(const index_t index, index_t dim) const {
        return complex->get(dim,index)->filtration();
	}

	void assemble_columns_to_reduce(std::vector<filtration_index_t>& simplices,
	                                std::vector<filtration_index_t>& columns_to_reduce,
	                                pivot_column_index_t& pivot_column_index, index_t dim, index_t num_simplices);

	void compute_dim_0_pairs(std::vector<filtration_index_t>& edges,
	                         std::vector<filtration_index_t>& columns_to_reduce) {
        //Get all edges and sort them
		filtered_union_find dset(complex->vertex_filtration());
		edges = get_edges();
		std::sort(edges.rbegin(), edges.rend(),
		          greater_filtration_or_smaller_index<filtration_index_t>());

		if(!python) { outfile << "# persistence intervals in dim 0:" << std::endl; }

		for (auto e : edges) {
			value_t birth = dset.link(complex->get(1,get_index(e))->vertices[0], complex->get(1,get_index(e))->vertices[1]);
			if (birth != -1) {
				if (get_filtration(e) > birth) {
					if(!python){
						outfile << " [" << birth << ", " << get_filtration(e) << ")" << std::endl;
					    outfile.flush();
					}
					else { finite_pairs[0].push_back(std::make_pair(birth,get_filtration(e))); }

				}
			} else {
				columns_to_reduce.push_back(e);
			}
		}
		std::reverse(columns_to_reduce.begin(), columns_to_reduce.end());

		for (index_t i = 0; i < n; ++i)
			if (dset.find(i) == i){
				if(!python){
					outfile << " [" << dset.filter(i) << ", )" << std::endl;
					outfile.flush();
				}
				infinite_pairs[0].push_back(dset.filter(i));
			}
	}

	template <typename Column, typename Iterator>
	filtration_entry_t add_coboundary_and_get_pivot(Iterator column_begin, Iterator column_end,
	                                              Column& working_coboundary, const index_t& dim,
                                                  std::priority_queue<filtration_entry_t, std::deque<filtration_entry_t>, smaller_index<filtration_entry_t>>&);

	void sort_columns(std::vector<filtration_index_t>& columns_to_reduce, index_t dimension) {
		std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
			greater_filtration_or_better_pivot_or_smaller_index(complex,dimension));
    }

	void compute_pairs(std::vector<filtration_index_t>& columns_to_reduce,
	                   pivot_column_index_t& pivot_column_index, index_t dim) {

		if(!python) { outfile << "# persistence intervals in dim " << dim << ":" << std::endl; }

        print_progress(dim,1);
        compressed_sparse_matrix<filtration_entry_t> reduction_coefficients;

		std::vector<filtration_entry_t> coface_entries;

		for (index_t index_column_to_reduce = 0; index_column_to_reduce < columns_to_reduce.size();
		     ++index_column_to_reduce) {
			auto column_to_reduce = columns_to_reduce[index_column_to_reduce];
            std::priority_queue<filtration_entry_t, std::deque<filtration_entry_t>, smaller_index<filtration_entry_t>> reduction_column;

			priority_queue_t<std::vector<filtration_entry_t>,
			                    greater_filtration_or_smaller_index<filtration_entry_t>>
			    working_coboundary(columns_to_reduce.size());

			value_t filtration = get_filtration(column_to_reduce);

#ifdef INDICATE_PROGRESS
			if ((index_column_to_reduce + 1) % 1 == 0)
				std::cout << "\033[K"
				          << "reducing column " << index_column_to_reduce + 1 << "/"
				          << columns_to_reduce.size() << " (filtration " << filtration << ")"
				          << std::flush << "\r";
#endif

			index_t index_column_to_add = index_column_to_reduce;

			filtration_entry_t pivot;

            reduction_coefficients.append_column();
			reduction_coefficients.push_back(filtration_entry_t(column_to_reduce, 1));

			while (true) {
                auto reduction_column_begin = reduction_coefficients.cbegin(index_column_to_add);
                auto reduction_column_end = reduction_coefficients.cend(index_column_to_add);

				pivot = add_coboundary_and_get_pivot(
				    reduction_column_begin, reduction_column_end,
				    working_coboundary, dim, reduction_column);

				if (get_index(pivot) > -1) {
                    auto pivot_column_idx = pivot_column_index[get_index(pivot)];

                    if (pivot_column_idx != INVALID_INDEX) {
                        index_column_to_add = pivot_column_idx;
                        continue;
                    } else {

						value_t death = get_filtration(pivot);
						if (death > filtration) {
							if(!python) {
								outfile << " [" << filtration << ", " << death << ")" << std::endl;
								outfile.flush();
							}
							else { finite_pairs[dim].push_back(std::make_pair(filtration,death)); }
						}

						pivot_column_index[get_index(pivot)] =  index_column_to_reduce;
                            reduction_coefficients.pop_back();
            				while (true) {
            					filtration_entry_t e = pop_pivot(reduction_column);
            					if (get_index(e) == -1) break;
            					reduction_coefficients.push_back(e);
            				}
						break;
					}
				} else if(get_index(pivot) == -1) {
					if(!python){
						outfile << " [" << filtration << ", )" << std::endl << std::flush;
						outfile.flush();
					}
					infinite_pairs[dim].push_back(filtration);
					break;
				}else {
					if(!python){
						outfile << "[?" << filtration << ", " << ", ?)" << std::endl;
						outfile.flush();
					}
					skipped_entries[dim]++;
					break;
				}
				// replace current column of reduction_coefficients (with a single diagonal 1 entry)
				// by reduction_column (possibly with a different entry on the diagonal)
			}
		}
#ifdef INDICATE_PROGRESS
		std::cout << "\033[K";
#endif
	}

	std::vector<filtration_index_t> get_edges();

	std::vector<value_t> num_infinite_pairs(){
		std::vector<value_t> out;
		for(auto i : infinite_pairs) out.push_back(i.size());
		return out;
	}

	void print_summary(){
		if(!python){
			std::vector<value_t> inf_pairs = num_infinite_pairs();
			outfile << std::endl;
			outfile << "# Betti Numbers:" << std::endl;
			for(index_t i = 0; i <= dim_max; i++){
				outfile << "#        dim H_" << i << " = " << inf_pairs[i];
				if( skipped_entries[i] > 0){
					outfile << " : (" << skipped_entries[i] << " entries skipped)";
				}
				outfile << std::endl;
			}
			outfile << std::endl;
			outfile << "# Cell Counts:" << std::endl;
			for(index_t i = 0; i <= dim_max; i++){
				outfile << "#        dim C_" << i << " = " << complex->number_of_cells(i) << std::endl;
			}
			outfile.flush();
			outfile.close();
		}
	}

	void compute_barcodes() {

        print_progress(0,2);
        print_progress(0,1);

		std::vector<filtration_index_t> simplices, columns_to_reduce;

		compute_dim_0_pairs(simplices, columns_to_reduce);

		for (index_t dim = 1; dim <= dim_max; ++dim) {
			pivot_column_index_t pivot_column_index(complex->number_of_cells(dim + 1), INVALID_INDEX);

			sort_columns(columns_to_reduce,dim);

			compute_pairs(columns_to_reduce, pivot_column_index, dim);

			if (dim < dim_max) {
				assemble_columns_to_reduce(simplices, columns_to_reduce, pivot_column_index,
				                           dim + 1, complex->number_of_cells(dim+1));
			}
		}
		print_summary();
	}


};

template <typename Column, typename Iterator>
filtration_entry_t tournser::add_coboundary_and_get_pivot(
    Iterator column_begin, Iterator column_end,
    Column& working_coboundary, const index_t& dim,
    std::priority_queue<filtration_entry_t, std::deque<filtration_entry_t>, smaller_index<filtration_entry_t>>& reduction_column) {
	index_t iterations = 0;
	for (auto it = column_begin; it != column_end; ++it) {
		filtration_entry_t simplex = *it;

        reduction_column.push(simplex);

		coface_entries.clear();
		simplex_coboundary_enumerator cofaces(simplex, dim, complex);
		while (cofaces.has_next()) {
			filtration_entry_t coface = cofaces.next();

            iterations++;
            working_coboundary.push(coface);
		}
		if (iterations > max_entries) {
			return filtration_entry_t(0,-2);
		}
	}

	return working_coboundary.get_pivot();
}

//returns a vector of all the edges where each is represented as a pair:
//(filtration,index), where the edge is the simplex at complex.get(1,index)
std::vector<filtration_index_t> tournser::get_edges() {
	std::vector<filtration_index_t> edges;
    int n = complex->number_of_cells(1);
	for ( index_t index = 0; index < n; index++) {
		edges.push_back(std::make_pair(complex->get(1,index)->filtration(), index));
	}
	return edges;
}

void tournser::assemble_columns_to_reduce(
    std::vector<filtration_index_t>& simplices, std::vector<filtration_index_t>& columns_to_reduce,
    pivot_column_index_t& pivot_column_index, index_t dim, index_t num_simplices ) {

	columns_to_reduce.clear();

	for (index_t index = 0; index < num_simplices; ++index) {
		if (pivot_column_index[index] == INVALID_INDEX) {
			value_t filtration = compute_filtration(index, dim);
			columns_to_reduce.push_back(std::make_pair(filtration, index));
		}
	}

	//std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
	//          greater_filtration_or_smaller_index<filtration_index_t>());
}


//END tournser
//-------------------------------------------------------------------------//
//BEGIN input parser

const arg_t parse_arguments(int argc, char** argv) {
	arg_t named_arguments;
	for (size_t i = 1; i < argc; ++i) {
		std::string arg(argv[i]);
		if (arg[0] == '-' && arg[1] == '-') {
			arg = arg.substr(2);
			named_arguments.insert(std::make_pair(arg, argv[++i]));
		}
	}
	return named_arguments;
}


//END input parser
//-------------------------------------------------------------------------//
//BEGIN main


int main(int argc, char** argv) {
	//input takes the form: ./tournser in_address out_address --approximate approx_val --filtration filt

    //read in arguments
    arg_t arguments = parse_arguments(argc,argv);
    arg_t::const_iterator it;
	const char* filename = argv[1];
	char* outname = argv[2];

    bool count_only = false;
    if ((it = arguments.find("count_only")) != arguments.end()){ count_only = true; }
    if (!count_only){
        std::ifstream f(outname);
    	if (f.good()){
    		std::cerr << "The output file already exists, aborting." << std::endl;
    		exit(-1);
    	} else { std::cout << "Printing output to : " << outname << std::endl; }
    }

	//initialise approximate functionality
	size_t max_entries = std::numeric_limits<size_t>::max();
	const char* filter = "vertexMax";
	if ((it = arguments.find("approximate")) != arguments.end()) {
		max_entries = atoi(it->second);
		std::cout << "Using Approximate Value: " << max_entries << std::endl;
	}
	if ((it = arguments.find("filtration")) != arguments.end()) { filter = it->second; }
	std::cout << "Using Filtration: " << filter << std::endl;

    bool dist = false;
    if ((it = arguments.find("print_dist")) != arguments.end()) { dist = true; }

	print_progress(0,4);

	delta_complex_t complex = read_graph(filename,filter,dist);
	complex.compute_oldest_cofaces();
	if ((it = arguments.find("print")) != arguments.end()){ complex.print_complex(it->second); }

    //create tournser object and compute persistent homology
    if (!count_only){
        tournser(&complex,outname,max_entries,false).compute_barcodes();
    }
}

//END main
//-------------------------------------------------------------------------//
