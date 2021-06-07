#include <iostream>
#include <vector>
#include<string>
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <bits/stdc++.h>
#include "spoa/include/spoa/spoa.hpp"
#include "alignment/alignment.h"
#include "bioparser/include/bioparser/fastq_parser.hpp"
#include "biosoup/include/biosoup/sequence.hpp"

using namespace std;

struct Sequence {  // or any other name
public:
	const char* name;
	std::uint32_t name_len;
      	const char* data;
	std::uint32_t data_len;
 //public:
  Sequence(  // required arguments
      const char* name, std::uint32_t name_len,
      const char* data, std::uint32_t data_len,
      const char* q, std::uint32_t q_len):
	name(name), name_len(name_len), data(data), data_len(data_len) {}
};

class Point{

private:
    int pointId, clusterId;
    int lenSequence;
    std::string seq;

public:
    Point(int id, std::string s){
        lenSequence = s.length();
        pointId = id;
        seq = s;
        clusterId = -1; //Initially not assigned to any cluster
    }

    int getlenSequence(){
        return lenSequence;
    }

    int getCluster(){
        return clusterId;
    }

    int getID(){
        return pointId;
    }

    void setCluster(int val){
        clusterId = val;
    }

    std::string getSeq() const{
        return this->seq;
    }
};

class Cluster{

private:
    int clusterId;
    Point centroid;
    std::vector<Point> points;

public:
	Cluster(int id, const Point& c):
		clusterId(id),
		centroid(c)
	{
		this->addPoint(c);
	}

    void addPoint(Point p){
        p.setCluster(this->clusterId);
        points.push_back(p);
    }

    void setCentroid(Point c){
        this->centroid = c;
    }

    std::string getCentroidValue(){
        std::string result = this->centroid.getSeq();
        return result;
    }

    bool removePoint(int pointId){
        int size = points.size();
        for(int i = 0; i < size; i++)
        {
            if(points[i].getID() == pointId)
            {
                points.erase(points.begin() + i);
                return true;
            }
        }
        return false;
    }

    int getId(){
        return clusterId;
    }

    Point getPoint(int pos){
        return points[pos];
    }

    int getSize(){
        return points.size();
    }

    Point getCentroid() {
        return centroid;
    }
};

std::string generateConsensus(std::vector<Point> &points){
	//std::cerr<<"generating consensus"<<std::endl;
	auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kOV, 1, 0, -1);  // linear gaps
	spoa::Graph graph{};
	for (const auto& it : points) {
		std::string el = it.getSeq();
		auto alignment = alignment_engine->Align(el, graph);
		graph.AddAlignment(alignment, el);
	}
	auto consensus = graph.GenerateConsensus();
	return consensus;
}

int alignPoints(Point p1, Cluster p2){
	std::string query = p1.getSeq();
	int query_len = query.length();
	std::string target = p2.getCentroidValue();
	int target_len = target.length();
	int minlen = -1;
	if(query_len < target_len){
		minlen = query_len;
	}
	else{
		minlen = target_len;
	}
	blonde::alignment::AlignmentType type = blonde::alignment::AlignmentType::kSemiGlobal;
	
	//std::string* cigar;
    	//unsigned int* target_begin;
	//orange::Alignment* a = new orange::Alignment(query, query_len, target, target_len, type, match, mismatch, gap);
    	//int result =  a->Align(query, query_len, target, target_len, type, match, mismatch, gap);
	int64_t result = blonde::alignment::Align(query.c_str(),query_len,target.c_str(),target_len,type, 1, 0, -1);
	//std::cerr<<"align result between "<<query<<" and "<<target<<": "<<result<<std::endl;
	//std::cerr<<std::endl;
	result = minlen - result;
	return result;
}

class KMeans{
private:
    int K, iters, total_points;
    //K = min similarity that allows point to get into cluster
    std::vector<Cluster> clusters;

    int getNearestClusterId(Point point, int K){
	//std::cout<<"in calc nearest cluster"<<std::endl;
        int NearestClusterId = -1;
        int temp = INT_MAX;
        int nClusters = clusters.size();
        for(int i = 0; i < nClusters; i++)
        {
            int dist = -1;
            dist = alignPoints(point, clusters[i]);
		//std::cout<<"calc dist: "<<dist<<std::endl;
            if(dist <= K && dist < temp) //manja udaljenost od K i od trenutno pronađenog
            {
		//std::cout<<"changed"<<std::endl;
                NearestClusterId = clusters[i].getId();
                temp = dist;
            }
        }

        return NearestClusterId;
    }

public:
    KMeans(int K, int iterations){
        this->K = K;
        this->iters = iterations;
    }

    void run(std::vector<Point>& all_points, int K){

        total_points = all_points.size();
	int tp = total_points;
	//generate starting cluster and align all sequences together
	auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kOV, 1, 0, -1);  // linear gaps
	spoa::Graph graph{};
	for (const auto& it : all_points) {
		std::string el = it.getSeq();
		auto alignment = alignment_engine->Align(el, graph);
		graph.AddAlignment(alignment, el);
	}
	auto consensus = graph.GenerateConsensus();
	std::cout<<"starting consensus: "<<consensus<<std::endl;
	total_points++;
	Point c = Point(total_points, consensus);
        Cluster cluster(0, c);
        clusters.push_back(cluster);
	std::vector<Point> tmp;
	auto msa = graph.GenerateMultipleSequenceAlignment();
	int i = 0;
	for (const auto& it : msa) {
		Point p = Point(i, it);
		i++;
		tmp.push_back(p);
	}
	all_points= tmp;
	int clusterCnt = 1;
        std::cout<<"Running K-Means Clustering.."<<std::endl;
        int iter = 1;
        while(true)
        {
            cout<<"Iter - "<<iter<<"/"<<iters<<endl;
            bool done = true;

            // Add all points to their nearest cluster
            for(int i = 0; i < tp; i++) 					//for all points
            {
		std::cout<<"point: "<<i<<std::endl;
                int currentClusterId = tmp[i].getCluster();
		//std::cout<<"current cluster id: "<< currentClusterId<<std::endl;
                int nearestClusterId = getNearestClusterId(tmp[i], K);
		//std::cout<<"nearest cluster id: "<< nearestClusterId<<std::endl;
                if(nearestClusterId == -1){
			//std::cout<<"generating new cluster"<<tmp[i].getSeq()<<std::endl;
                	Cluster cluster(clusterCnt, tmp[i]);
			clusterCnt++;
                	clusters.push_back(cluster);
                }
                else{
                    if(currentClusterId != nearestClusterId)
                    {
			//std::cout<<"checkpoint 1"<<std::endl;
                        if(currentClusterId != -1){
                            for(int j=0; j< clusters.size(); j++){
                            	if(clusters[j].getId() == currentClusterId){
                                    	clusters[j].removePoint(tmp[i].getID());//if it is already in other cluster, remove it
					//std::cout<<"removed"<<std::endl;
                                }
                            }
                        }

                        for(int j=0; j< clusters.size(); j++){
                            if(clusters[j].getId() == nearestClusterId){
				//add point to new cluster
                                clusters[j].addPoint(tmp[i]);
				//std::cout<<"added"<<std::endl;
				//recalculate consensus
				int ClusterSize = clusters[j].getSize();
				vector<Point> newPoints;
				for(int k = 0; k < ClusterSize; k++){
				    newPoints.push_back(clusters[j].getPoint(k));
				}
				std::string consensus = generateConsensus(newPoints);
				//std::cout<<"new consensus: "<<consensus<<std::endl;
				int tmpCnt = total_points + clusterCnt;
				Point c = Point(tmpCnt, consensus);
				clusters[j].setCentroid(c);
				break;
		            }
                        }
			tmp[i].setCluster(nearestClusterId);
		        done = false;
                   }
		}
            }
		/*
            // Recalculating centroid of each cluster
		for(int i = 0; i < clusters.size(); i++)
		{
		    int ClusterSize = clusters[i].getSize();
		    if(ClusterSize > 0){
		        vector<Point> newPoints;
		        for(int j = 0; j < ClusterSize; j++){
		            newPoints.push_back(clusters[i].getPoint(j));
		        }
		        //calculate consensus
		        std::string consensus = generateConsensus(newPoints);
		        total_points ++;
		        Point c = Point(total_points, consensus);
		        clusters[i].setCentroid(c);
		    }
		}*/
            if(done || iter >= iters)
            {
                std::cout << "Clustering completed in iteration : " <<iter<<std::endl<<std::endl;
                break;
            }
            iter++;
        }
        cout<<"========================"<<endl<<endl;

        //Write cluster centers to file
        ofstream outfile;
        outfile.open("clusters.txt");
	std::cerr<<"Number of clusters:"<<clusters.size()<<std::endl;
        if(outfile.is_open()){
            for(int i=0; i<clusters.size(); i++){
                std::cout<<"Cluster "<<i<<" centroid : "<<clusters[i].getCentroidValue()<<std::endl;
		std::cout<<"Number of points in this cluster: "<<clusters[i].getSize()<<std::endl;
                outfile<<"Cluster "<<i<<" centroid : "<<clusters[i].getCentroidValue()<<std::endl;
		outfile<<"Number of points in this cluster: "<<clusters[i].getSize()<<std::endl;
                std::cout<<std::endl;
                outfile<<std::endl;
            }
            outfile.close();
        }
        else{
            cout<<"Error: Unable to write to clusters.txt";
        }

    }
};
/*
std::vector<std::string> loadFastq(std::string fileName){
	std::vector<std::string> result;
	auto fragment_parser = bioparser::Parser<Sequence>::Create<bioparser::FastqParser>(fileName);
        // parse in chunks
        std::vector<std::unique_ptr<Sequence>> fragments;
        std::uint32_t chunk_size = 500 * 1024 * 1024;  // 500 MB
        for (auto t = fragment_parser->Parse(chunk_size); !t.empty(); t = fragment_parser->Parse(chunk_size)) {
            fragments.insert(fragments.end(), std::make_move_iterator(t.begin()), std::make_move_iterator(t.end()));
        }
	for(int i = 0; i < fragments.size();i++){
		std::string el = ((char*)((fragments[i])->data)).c_str();
		result.push_back(el);
	}
	return result;

}*/
bool getFileContent(std::string fileName, std::vector<std::string> & vecOfStrs)
{
    std::ifstream in(fileName.c_str());
    if(!in)
    {
        std::cerr << "Cannot open the File : "<<fileName<<std::endl;
        return false;
    }
    std::string str;
    while (std::getline(in, str))
    {
        if(str.size() > 0)
            vecOfStrs.push_back(str);
    }
    in.close();
    return true;
}



int main(int argc, char** argv) {
	
	std::vector<std::string> sequences_from_file;
	bool result = getFileContent("seq_j29B.txt", sequences_from_file);
	
	std::vector<Point> points;
	int i = 0;
	//get all points of size 296
	for(const auto& it : sequences_from_file){
		if(it.length() > 293 && it.length() < 299){
			Point p = Point(i, it);
			points.push_back(p);
		}
		i++;
	}
	int K = 6; //granicna duljina
	int iters = 100;
	std::cout<<points.size();
	KMeans kmeans(K, iters);
	kmeans.run(points, K);

	return 0;
}
