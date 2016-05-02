/*
	Dijkstra and A* C++ implementation
	CMPT 307, Summer 2015
	Byron Chiang 2*******4
	Braeden Mulligan 3*******9
	
	See the readme file for more information
*/

#include <queue>    // STL class priority queue library
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <string>
#include <iomanip>  // for increasing display precision of doubles
#include <limits>   // for setting infinity values
#include <utility>  // for creating pair/tuple
#include <cmath>
#include <ctime>
#include <algorithm>    // max function only... did not call some crazy algorithm
using namespace std;

ofstream output;

// define graph G
// index of is is the node #
typedef pair<int, double> neighbor_cost;  // adjacent neighbor node id (int) and cost (double) to get there
typedef vector<neighbor_cost> set_neighbors; // set of adjacent neighbors linked to a node
typedef pair <double, double> coordinates;
typedef struct vertex_t // create a vertex type
{
    coordinates lat_long; // coordinates of the node
    set_neighbors adj_id; // a vector to store adjacent list of neighbors and costs to get there
} Vertex;

typedef vector<Vertex> graph;    // graph containing a set of the above

// Record nodes visited by each algorithm for comparison.
int sum_dijkstra = 0;
int sum_astar = 0;
int sum_astar_landmark = 0;

// comparator to make priority queue pop minimum element
struct costCompare
{
    bool operator()(const neighbor_cost &c1, const neighbor_cost &c2) const
    {
        return c1.second > c2.second;
    }
};


// calculates the distance between two sets of coordinates
// param:   lat_long_1: latitude and longitude of vertex 1
//          lat_long_2: latitude and longitude of vertex 2
double coordDist( coordinates lat_long_1, coordinates lat_long_2 )
{
    // otherwise
    double constant = M_PI / 180;   // convert degrees to radians (Binay's eq = 2pi/360)
    double dlat = constant * (lat_long_2.first - lat_long_1.first);
    double mlat = constant * (lat_long_2.first + lat_long_1.first)/2;
    double dlon = constant * (lat_long_2.second - lat_long_1.second);
    return 6371009 * pow(pow(dlat, 2.0) + pow(cos(mlat)*dlon, 2.0), .5);
}


// input data from input file
// param:   G: graph
void inputData(graph &G)
{
    int linecount = 1;
    string str = "";
    istringstream iss;

    // input data from "graph1000.txt"
    ifstream input;
    input.open("graph1000.txt");
    if (!input)
    {
        cerr << " failed to read graph1000.txt" << endl;
        exit(1);
    }


    // get the vertex id and its coordinates
    int i = 0;
    while(linecount <= G.size())
    {
        input.ignore(10, ':');
        input.ignore(1, ' ');

        getline(input, str, ',');
        double latitude = atof(str.c_str());
        input.ignore(1);

        getline(input, str);
        double longitude = atof(str.c_str());

        G[i].lat_long = make_pair(latitude, longitude);

//        cout << "Just imported node #" << i;
//        cout << " with x = " << setprecision(15) << G[i].lat_long.first << " and ";
//        cout << "y = " << setprecision(15) << G[i].lat_long.second << endl;

        linecount++;
        i++;
    }


    // go back and get the adjacency list
    i = 0;
    input.ignore(1, '\n');  // ignore line 1001
    while (i < 1000)
//    while (i < 10)
    {
        input.ignore(10, ':');
        input.ignore(1, ' ');


        getline(input, str);
        iss.str(str);
        int adj_node;
        while (iss >> adj_node)
        {
            if (iss.peek() == ',')
            iss.ignore();

            // insert with adjusted node position to accomodate shifted vertex index by -1
            double length = coordDist(G[i].lat_long, G[adj_node - 1].lat_long);
            G[i].adj_id.push_back( make_pair(adj_node - 1, length) );
        }

        iss.clear();

        i++;
    }

//    for (int i = 0; i < G.size(); i++)
//    {
//        cout << i << " has " << G[i].adj_id.size() << " adjacent nodes" << endl;
//        for (int j = 0; j < G[i].adj_id.size(); j++)
//            cout << i << " is adjacent to " << G[i].adj_id[j].first << endl;
//    }

    input.close();  // close file
}


// output graph just imported - confirms that data is accurately imported
// param:   G: graph
void outputData(graph &G)
{
    ofstream output;
    output.open("output_graph.txt");
    if (!output)
    {
        cerr << "failed to write to output.txt" << endl;
        exit(1);
    }

    output << "The current size of G is " << G.size() << endl;
    for (int i = 0; i < G.size(); i++)
    {
		output << "number of adjacent vertices at node " << i << " is " << G[i].adj_id.size() << endl;
        for (int j = 0; j < G[i].adj_id.size(); j++)
        {
            output << "it is adjacent to node " << G[i].adj_id[j].first << " which has a weight of " << setprecision(15) << G[i].adj_id[j].second << endl;
        }
		output << endl;
    }

    output << "-------------------------------- OUTPUT FINISHED --------------------------------";

    output.close(); // close file
}


// A* algorithm as per pseudocode
// pre: G exists with positive edges
// param:   G - graph
//          s - node id to start A*
//          t - node id to end A*
// post: returns the shortest path
double astar(const graph &G, const int s, const int t)
{
    int counter = 0;

    // set up ofstream.  output data to "output_astar.txt"
    ofstream output;
    output.open("output_astar.txt", ios::app);
    if (!output)
    {
        cerr << "failed to write to output_astar.txt" << endl;
        exit(1);
    }


    vector<double> dist(G.size());
    vector<double> est(G.size());
    vector<int> pred(G.size());
    vector<int> visited;


    output << "------------------------------------- A* BEGINS -------------------------------------" << endl;
    output << "starting node s = " << s << endl;
    output << "destination node t = " << t << endl;

    for (int u = 0; u < G.size(); u++)
    {
        dist[u] = numeric_limits<double>::infinity(); // set all distances to infinity
        est[u] = numeric_limits<double>::infinity(); // set all estimates to infinity
        pred[u] = -1;   // came from predecessor node.  -1 = nil
    }


    dist[s] = (double)0;
    double rem_s = coordDist(G[s].lat_long, G[t].lat_long);
    est[s] = dist[s] + rem_s;
    output << "the initial estimate from s to t is " << setprecision(15) << est[s] << endl;
    output << "Please note that all nodes have been offset by 1 (node range is not 1-1000 but 0-999)" << endl;
    output << "-------------------------------------------------------------------------------------" << endl;


    // initialize a priority queue that sets minimum element = highest priority
    priority_queue<neighbor_cost, vector<neighbor_cost>, costCompare > H;

    H.push(make_pair(s, est[s]));

    while (!H.empty())
    {
        // set u = extractMin
        int u = H.top().first; // the shortest neighbor in the queue
        H.pop();

        if (u == t)
        {
            output << "-------------------------------- A*'S RESULTS --------------------------------" << endl;
            output << "the total path distance from s: " << setprecision(15) << dist[u] << " meters" << endl;
            output << "nodes visited: " << counter << endl;

			sum_astar = sum_astar + counter;
            output << "------------------------------------------------------------------------------" << endl << endl;
            output.close(); // close file

            return dist[u];
        }


        vector<int>::iterator found;
//        cout << "current size of visited is " << visited.size() << endl;
        found = find(visited.begin(), visited.end(), u);    // find if v has aleady been visited
        if (found == visited.end()) // if item has not been visited
        {
            for (int i = 0; i < G[u].adj_id.size(); i++)   // for all edges e of E
            {
                int v = G[u].adj_id[i].first;  // neighbor node id
                double c = G[u].adj_id[i].second;  // distance from u to its neighbor

                // relax portion
                if (dist[v] > dist[u] + c)
                {
                    dist[v] = dist[u] + c;
                    pred[v] = u;
                    double rem_v = coordDist(G[v].lat_long, G[t].lat_long);
                    est[v] = dist[v] + rem_v;
                    H.push(make_pair(v, est[v]));   // decreaseKey
                }
            }

            visited.push_back(u);
//            output << "visited " << u << endl;
            counter++;
        }

        // do nothing if already visited
    }


    output << "-------------------------------- A*'S RESULTS --------------------------------" << endl;
    output << "t is unreachable.  dist[t] = " << dist[t] << endl;
    output << "nodes visited: " << counter << endl;
    output << "------------------------------------------------------------------------------" << endl << endl;
    output.close(); // close file

    return dist[t];
}



// generate three landmarks in random
// pre: size of random_node has to be three
//      - because I am too lazy to write a comparator for this.  Used max function instead
// param:   G - graph containing all the imported data
//          startn - start point (node#)
//          endn - end point (node#)
//          random_node - vector containing randomly generated landmarks (node #s)
double landmark(const graph &G,
                int startp,
                int endp,
                const vector<int> &random_nodes)
{
    vector<neighbor_cost> rem_lm;
    vector<neighbor_cost> dist_lm_startp;
    vector<neighbor_cost> dist_lm_endp;


    // initialize shortest path distances from s to landmarks
    for (int i = 0; i < 3; i++)
    {
        int node = random_nodes[i];
        dist_lm_startp.push_back( make_pair(node, coordDist(G[node].lat_long, G[startp].lat_long)) );
        dist_lm_endp.push_back( make_pair(node, coordDist(G[node].lat_long, G[endp].lat_long)) );
        rem_lm.push_back( make_pair(node, abs(dist_lm_startp[i].second - dist_lm_endp[i].second)) );
    }

    return max( max(rem_lm[0].second, rem_lm[1].second), rem_lm[2].second);

}



// A* algorithm with landmark
// pre: G exists with positive edges
// param:   G - graph
//          s - node id to start A*
//          t - node id to end A*
// post: returns the shortest path
double astar_landmark(const graph &G, const int s, const int t)
{
    int counter = 0;

    // set up ofstream.  output data to "output_astar_landmark.txt"
    ofstream output;
    output.open("output_astar_landmark.txt", ios::app);
    if (!output)
    {
        cerr << "failed to write to output_astar_landmark.txt" << endl;
        exit(1);
    }


    vector<double> dist(G.size());
    vector<double> est(G.size());
    vector<int> pred(G.size());
    vector<int>random_nodes(3);
    vector<int> visited;


    output << "--------------------------------- A* LANDMARK BEGINS --------------------------------" << endl;
    output << "starting node s = " << s << endl;
    output << "destination node t = " << t << endl;

    for (int u = 0; u < G.size(); u++)
    {
        dist[u] = numeric_limits<double>::infinity(); // set all distances to infinity
        est[u] = numeric_limits<double>::infinity(); // set all estimates to infinity
        pred[u] = -1;   // came from predecessor node.  -1 = nil
    }

    dist[s] = (double)0;

    // generate three random nodes
    for (int i = 0; i < random_nodes.size(); i++)
        random_nodes[i] = rand() % 1000;

    // apply landmark estimate
    est[s] = dist[s] + landmark(G, s, t, random_nodes);

    output << "the initial estimate from s to t is " << setprecision(15) << est[s] << endl;
    output << "Please note that all nodes have been offset by 1 (node range is not 1-1000 but 0-999)" << endl;
    output << "-------------------------------------------------------------------------------------" << endl;



    // initialize a priority queue that sets minimum element = highest priority
    priority_queue<neighbor_cost, vector<neighbor_cost>, costCompare > H;

    H.push(make_pair(s, est[s]));

    while (!H.empty())
    {
        // set u = extractMin
        int u = H.top().first; // the shortest neighbor in the queue
        H.pop();

        if (u == t)
        {
            output << "-------------------------------- A* LANDMARK'S RESULTS --------------------------------" << endl;
            output << "the total path distance from s: " << setprecision(15) << dist[u] << " meters" << endl;
            output << "nodes visited: " << counter << endl;

			sum_astar_landmark = sum_astar_landmark + counter;
            output << "---------------------------------------------------------------------------------------" << endl << endl;
            output.close(); // close file

            return dist[u];
        }



        vector<int>::iterator found;
//        cout << "current size of visited is " << visited.size() << endl;
        found = find(visited.begin(), visited.end(), u);    // find if v has aleady been visited
        if (found == visited.end()) // if item has not been visited
        {
            for (int i = 0; i < G[u].adj_id.size(); i++)   // for all edges e of E
            {
                int v = G[u].adj_id[i].first;  // neighbor node id
                double c = G[u].adj_id[i].second;  // distance from u to its neighbor

                // relax portion
                if (dist[v] > dist[u] + c)
                {
                    dist[v] = dist[u] + c;
                    pred[v] = u;
                    est[v] = dist[v] + landmark(G, v, t, random_nodes);
                    H.push(make_pair(v, est[v]));   // decreaseKey
                }
            }

            visited.push_back(u);
//            output << "visited " << u << endl;
            counter++;
        }

        // do nothing if already visited
    }


    output << "-------------------------------- A* LANDMARK'S RESULTS --------------------------------" << endl;
    output << "t is unreachable.  dist[t] = " << dist[t] << endl;
    output << "nodes visited: " << counter << endl;
    output << "---------------------------------------------------------------------------------------" << endl << endl;
    output.close(); // close file

    return dist[t];
}



// Dijkstra's shortest path algorithm as per pseudocode
// pre: G exists with positive edges
// param:   G - graph
//          s - node id to start dijkstra
//          t - node id to end dijkstra
// post: returns the shortest path
double dijkstra(const graph &G, const int s, const int t)
{
    int counter = 0;
    // set up ofstream.  output data to "output_dijkstra.txt"
    ofstream output;
    output.open("output_dijkstra.txt", ios::app);
    if (!output)
    {
        cerr << "failed to write to output_dijkstra.txt" << endl;
        exit(1);
    }


    // Dijkstra stuff
    vector<double> dist(G.size());
    vector<int> pred(G.size());
    vector<int> visited;

    output << "---------------------------------- DIJKSTRA BEGINS ----------------------------------" << endl;
    output << "starting node s = " << s << endl;
    output << "destination node t = " << t << endl;

    for (int u = 0; u < G.size(); u++)
    {
        dist[u] = numeric_limits<double>::infinity(); // set all distances to infinity
        pred[u] = -1;   // came from predecessor node.  -1 = nil
    }

    dist[s] = (double)0;
    output << "Please note that all nodes have been offset by 1 (node range is not 1-1000 but 0-999)" << endl;

    // initialize a priority queue that sets minimum element = highest priority
    priority_queue<neighbor_cost, vector<neighbor_cost>, costCompare > H;

    H.push(make_pair(s, dist[s]));
    output << "node #" << s << " and " << "a starting distance " << setprecision(15) << dist[s] << " pushed into the queue" << endl;

    while (!H.empty())
    {
        // set u = extractMin
        int u = H.top().first; // the shortest neighbor in the queue
        H.pop();

        if (u == t)
        {

            output << "-------------------------------- DIJKSTRA'S RESULTS --------------------------------" << endl;
            output << "reached final destination t: " << u << endl;
            output << "the total path length from s to t: " << setprecision(15) << dist[u] << " meters" << endl;
            output << "nodes visited: " << counter << endl;

			sum_dijkstra = sum_dijkstra + counter;
            output << "------------------------------------------------------------------------------------" << endl << endl;
            output.close(); // close file
            return dist[u];
        }


        vector<int>::iterator found;
//        cout << "current size of visited is " << visited.size() << endl;
        found = find(visited.begin(), visited.end(), u);    // find if v has aleady been visited
        if (found == visited.end()) // if item has not been visited
        {
            for (int i = 0; i < G[u].adj_id.size(); i++)   // for all edges e of E
            {
                int v = G[u].adj_id[i].first;  // neighbor node id
                double c = G[u].adj_id[i].second;  // distance from u to its neighbor

                // relax portion
                if (dist[v] > dist[u] + c)
                {
                    dist[v] = dist[u] + c;
                    pred[v] = u;

                    H.push(make_pair(v, dist[v]));   // decreaseKey
                }
            }

            visited.push_back(u);
//            output << "visited " << u << endl;
            counter++;
        }

        // do nothing if already visited
    }


    output << "-------------------------------- DIJKSTRA COMPLETES --------------------------------" << endl;
    output << "t is unreachable.  dist[t] = " << dist[t] << endl;
    output << "nodes visited: " << counter << endl;
    output << "------------------------------------------------------------------------------------" << endl;
    output.close(); // close file


    return dist[t];
}




int main()
{
	// Clear old outputs.
	remove("output_dijkstra.txt");
	remove("output_astar.txt");
	remove("output_astar_landmark.txt");

    // initialize a Graph containing all 1000 vertexes
    int numPoints = 1000;
    graph G(numPoints);
    inputData(G);

    outputData(G);

//    double result = dijkstra(G, 1, 10);
//    cout << "the result is: " << setprecision(15) << result << endl;
//
//    double result2 = astar(G, 1, 10);
//    cout << "the result is: " << setprecision(15) << result2 << endl;
//
//    double result3 = astar_landmark(G, 1, 10);
//    cout << "the result is: " << setprecision(15) << result2 << endl;

    // initialize random number generator
    srand(time(0));

    // generate 20 random queries and apply Dijkstra, A*, and A*  Landmark
    for (int i = 0; i < 20; i++)
    {
        int random_s = rand() % 1000;
        int random_t = rand() % 1000;
        dijkstra(G, random_s, random_t);
		
        astar(G, random_s, random_t);
        astar_landmark(G, random_s, random_t);
    }

	float average_dijkstra = (float)sum_dijkstra / 20.0;
	float average_astar = (float)sum_astar / 20.0;
	float average_astar_landmark = (float)sum_astar_landmark / 20.0;

	ofstream output;
	output.open("output_comparison.txt");

	output << "Dijkstra visited an average of " << average_dijkstra << " nodes." << endl;
	output << "A* visited an average of " << average_astar << " nodes." << endl;
	output << "A* with landmarks visited an average of " << average_astar_landmark << " nodes." << endl << endl;

	output << ("Please see the other output files for the shortest paths and number of nodes visited after each randomized run of the algorithms.") << endl;
	
    return 0;
}

