#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
// #include <bits/stdc++.h>
#include "spoa/include/spoa/spoa.hpp"
#include "alignment/alignment.h"
#include "bioparser/include/bioparser/fastq_parser.hpp"
#include "biosoup/include/biosoup/sequence.hpp"

using namespace std;

struct Sequence
{ 
public:
    const char *name;
    std::uint32_t name_len;
    const char *data;
    std::uint32_t data_len;
    
    Sequence( 
        const char *name, std::uint32_t name_len,
        const char *data, std::uint32_t data_len,
        const char *q, std::uint32_t q_len) : name(name), name_len(name_len), data(data), data_len(data_len)
    {
    }
};

class Point
{

private:
    int pointId, clusterId;
    int lenSequence;
    std::string seq;

public:
    Point(int id, std::string s)
    {
        lenSequence = s.length();
        pointId = id;
        seq = s;
        clusterId = -1; //Initially not assigned to any cluster
    }

    int getlenSequence()
    {
        return lenSequence;
    }

    int getCluster()
    {
        return clusterId;
    }

    int getID()
    {
        return pointId;
    }

    void setCluster(int val)
    {
        clusterId = val;
    }

    std::string getSeq() const
    {
        return this->seq;
    }
};

class Cluster
{

private:
    int clusterId;
    Point centroid;
    std::vector<Point> points;

public:
    Cluster(int id, const Point &c) : clusterId(id),
                                      centroid(c)
    {
        this->addPoint(c);
    }

    void addPoint(Point p)
    {
        p.setCluster(this->clusterId);
        points.push_back(p);
    }

    void setCentroid(Point c)
    {
        this->centroid = c;
    }

    std::string getCentroidValue()
    {
        std::string result = this->centroid.getSeq();
        return result;
    }

    bool removePoint(int pointId)
    {
        int size = points.size();
        for (int i = 0; i < size; i++)
        {
            if (points[i].getID() == pointId)
            {
                points.erase(points.begin() + i);
                return true;
            }
        }
        return false;
    }

    int getId()
    {
        return clusterId;
    }

    Point getPoint(int pos)
    {
        return points[pos];
    }

    int getSize()
    {
        return points.size();
    }

    Point getCentroid()
    {
        return centroid;
    }
};

std::string generateConsensus(std::vector<Point> &points)
{
    auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kOV, 1, 0, -1);
    spoa::Graph graph{};
    for (const auto &it : points)
    {
        std::string el = it.getSeq();
        auto alignment = alignment_engine->Align(el, graph);
        graph.AddAlignment(alignment, el);
    }
    auto consensus = graph.GenerateConsensus();
    return consensus;
}

int alignPoints(Point p1, Cluster p2)
{
    std::string query = p1.getSeq();
    int query_len = query.length();
    std::string target = p2.getCentroidValue();
    int target_len = target.length();
    int minlen = -1;
    if (query_len < target_len)
    {
        minlen = query_len;
    }
    else
    {
        minlen = target_len;
    }
    alignment::AlignmentType type = alignment::AlignmentType::kSemiGlobal;

    int64_t result = alignment::Align(query.c_str(), query_len, target.c_str(), target_len, type, 1, 0, -1);
    result = minlen - result;
    return result;
}

class KMeans
{
private:
    int K, iters, total_points;
    //K = min similarity that allows point to get assigned to cluster
    std::vector<Cluster> clusters;

    int getNearestClusterId(Point point, int K)
    {
        int NearestClusterId = -1;
        int temp = INT_MAX;
        int nClusters = clusters.size();
        for (int i = 0; i < nClusters; i++)
        {
            int dist = -1;
            dist = alignPoints(point, clusters[i]);
            if (dist <= K && dist < temp) //lower distance than K and current one
            {
                NearestClusterId = clusters[i].getId();
                temp = dist;
            }
        }
        return NearestClusterId;
    }

public:
    KMeans(int K, int iterations)
    {
        this->K = K;
        this->iters = iterations;
    }

    void run(std::vector<Point> &all_points, int K)
    {

        total_points = all_points.size();
        int tp = total_points;
        //generate starting cluster and align all sequences together
        auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kOV, 1, 0, -1);
        spoa::Graph graph{};
        for (const auto &it : all_points)
        {
            std::string el = it.getSeq();
            auto alignment = alignment_engine->Align(el, graph);
            graph.AddAlignment(alignment, el);
        }
        auto consensus = graph.GenerateConsensus();
        total_points++;
        Point c = Point(total_points, consensus);
        Cluster cluster(0, c);
        clusters.push_back(cluster);
        std::vector<Point> tmp;
        auto msa = graph.GenerateMultipleSequenceAlignment();
        int i = 0;
        for (const auto &it : msa)
        {
            Point p = Point(i, it);
            i++;
            tmp.push_back(p);
        }
        all_points = tmp;
        int clusterCnt = 1;
        std::cout << "Starting algorithm for K = "<< K << std::endl;
        int iter = 1;
        while (true)
        {
            cout << "Iter - " << iter << "/" << iters << "             " << endl;
            bool done = true;

            // Add all points to their nearest cluster
            for (int i = 0; i < tp; i++) //for all points
            {
                std::cout << "Point: " << i << "/" << tp << "\r";
                std::cout.flush();
                int currentClusterId = tmp[i].getCluster();
                int nearestClusterId = getNearestClusterId(tmp[i], K);
                if (nearestClusterId == -1)
                {
                    Cluster cluster(clusterCnt, tmp[i]);
                    clusterCnt++;
                    clusters.push_back(cluster);
                }
                else
                {
                    if (currentClusterId != nearestClusterId)
                    {
                        if (currentClusterId != -1)
                        {
                            for (int j = 0; j < clusters.size(); j++)
                            {
                                if (clusters[j].getId() == currentClusterId)
                                {
                                    clusters[j].removePoint(tmp[i].getID()); //if it is already in other cluster, remove it
                                }
                            }
                        }

                        for (int j = 0; j < clusters.size(); j++)
                        {
                            if (clusters[j].getId() == nearestClusterId)
                            {
                                //add point to new cluster
                                clusters[j].addPoint(tmp[i]);
                                //recalculate consensus
                                int ClusterSize = clusters[j].getSize();
                                vector<Point> newPoints;
                                for (int k = 0; k < ClusterSize; k++)
                                {
                                    newPoints.push_back(clusters[j].getPoint(k));
                                }
                                std::string consensus = generateConsensus(newPoints);
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
            if (done || iter >= iters)
            {
                std::cout << "Clustering completed in iteration : " << iter << std::endl
                          << std::endl;
                break;
            }
            iter++;
        }
        std::cout << "Finished." << std::endl << std::endl;

        //Write cluster centers to file
        ofstream outfile;
        std::string fileName;
        fileName = "clusters_" + std::to_string(K) + ".txt";
        outfile.open(fileName);
        std::cerr << "Number of clusters:" << clusters.size() << std::endl;
        if (outfile.is_open())
        {
            outfile << "Performing for K = " + std::to_string(K) << std::endl
                    << std::endl;
            for (int i = 0; i < clusters.size(); i++)
            {
                std::cout << "Cluster " << i << " centroid : " << clusters[i].getCentroidValue() << std::endl;
                std::cout << "Number of points in this cluster: " << clusters[i].getSize() << std::endl;
                outfile << "Cluster " << i << " centroid : " << clusters[i].getCentroidValue() << std::endl;
                outfile << "Number of points in this cluster: " << clusters[i].getSize() << std::endl;
                std::cout << std::endl;
                outfile << std::endl;
            }
            outfile.close();
        }
        else
        {
            cout << "Error: Unable to write to clusters.txt";
        }
    }
};

bool getFileContent(std::string fileName, std::vector<std::string> &vecOfStrs)
// reads vectors from the file
{
    std::ifstream in(fileName.c_str());
    if (!in)
    {
        std::cerr << "Cannot open the File : " << fileName << std::endl;
        return false;
    }
    std::string str;
    while (std::getline(in, str))
    {
        if (str.size() > 0)
            vecOfStrs.push_back(str);
    }
    in.close();
    return true;
}

void runIteration(int K, int iters, std::vector<Point> &points)
// runs single iteration, used to dispose objects after one iteration
{
    KMeans kmeans(K, iters);
    kmeans.run(points, K);
}

int main(int argc, char **argv)
{

    std::vector<std::string> sequences_from_file;
    //currently running for file J30_B_CE_IonXpress.fastq
    bool result = getFileContent("seq_j30B.txt", sequences_from_file);

    std::vector<Point> points;
    int i = 0;
    //get all points of size 296 +/- 3
    for (const auto &it : sequences_from_file)
    {
        if (it.length() > 293 && it.length() < 299)
        {
            Point p = Point(i, it);
            points.push_back(p);
        }
        i++;
    }
    //int K = 17; //limiting distance
    int iters = 100;
    for (int j = 20; j > 15; j--)
    {
        runIteration(j, iters, points);
    }

    return 0;
}
