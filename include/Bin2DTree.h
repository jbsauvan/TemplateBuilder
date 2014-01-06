

#ifndef BIN2DTREE_H
#define BIN2DTREE_H

#include <iostream>
#include <vector>
#include <algorithm>
#include "TLine.h"
#include "TH2F.h"

class Bin2DLeaf
{
    /* Leaf of a Bin2DTree.
    It encodes a 2D bin and stores the entries contained in this bin: (x,y) coordinates + event weight.
    Each leaf is identified with an index.
    */
    public:
        Bin2DLeaf();
        Bin2DLeaf(double xmin, double xmax, double ymin, double ymax);
        ~Bin2DLeaf();

        bool operator <(const Bin2DLeaf& leaf) const
        {
            return m_index < leaf.index();
        }

        void setBinBoundaries(double xmin, double xmax, double ymin, double ymax);
        const std::vector< std::pair<double,double> >& getBinBoundaries();
        double getXMin();
        double getXMax();
        double getYMin();
        double getYMax();
        double getXWidth();
        double getYWidth();
        double getXCenter();
        double getYCenter();
        bool isNeighbor(Bin2DLeaf* leaf);
        unsigned int getNEntries();
        double getSumOfWeights();
        const std::vector< std::vector<double> >& getEntries();
        std::vector<double> percentiles(const std::vector<double>& q, unsigned int axis=0);
        double densityGradient(unsigned int axis=0, double q=10.);
        double density(double xmin, double xmax, unsigned int axis=0);
        bool inBin(double xi, double yi);
        bool addEntry(double xi, double yi, double wi);
        std::vector<TLine*> getBoundaryTLines();

        unsigned int index() const {return m_index;}
        void setIndex(unsigned int index) {m_index = index;}

    private:
        unsigned int m_index;
        std::vector< std::pair<double,double> > m_binBoundaries;
        std::vector< std::vector<double> > m_entries;


};


class Bin2DTree
{
    public:
        Bin2DTree(double xmin, double xmax, double ymin, double ymax, const std::vector< std::vector<double> >& entries);
        ~Bin2DTree();
        void addEntry(double xi, double yi, double wi);
        std::vector< std::pair<double,double> > getBinBoundaries();
        double getXMin();
        double getXMax();
        double getYMin();
        double getYMax();
        unsigned int getNEntries();
        double getSumOfWeights();
        std::vector< std::vector<double> > getEntries();
        Bin2DLeaf* getLeaf(double xi, double yi);
        std::vector<Bin2DLeaf*> getLeaves();
        std::vector<Bin2DTree*> getTerminalNodes();
        std::vector<Bin2DLeaf*> findNeighborLeaves(Bin2DLeaf* leaf);
        unsigned int getNLeaves();
        unsigned int maxLeafIndex();
        unsigned int minLeafEntries() const {return m_minLeafEntries;}
        double maxAxisAsymmetry() const {return m_maxAxisAsymmetry;}
        void build();
        std::vector<TLine*> getBoundaryTLines();
        TH2F* fillHistogram();
        std::pair<TH2F*,TH2F*> fillWidths();

        Bin2DLeaf* leaf(){return m_leaf;}
        void setGridConstraint(TH2F* gridConstraint);
        void setVetoSplitXY(unsigned int axis, bool veto){m_vetoSplitXY[axis] = veto;}
        void setMinLeafEntries(unsigned int minLeafEntries){m_minLeafEntries = minLeafEntries;}
        void setMaxAxisAsymmetry(double maxAxisAsymmetry){m_maxAxisAsymmetry = maxAxisAsymmetry;}
        bool vetoSplitXY(unsigned int axis){return m_vetoSplitXY[axis];}

    private:
        void splitLeaf(double cut, unsigned int maxLeafIndex, unsigned int axis=0);
        void findBestSplit(Bin2DTree*& bestNode, unsigned int& axis, double& gradient);
        void constrainSplit(int axis, double& cut, bool& veto);
        void minimizeLongBins(Bin2DTree* tree, unsigned int axis, double& cut, bool& veto);

        std::vector<Bin2DTree*> m_treeSons;
        unsigned int m_cutAxis;
        double m_cut;
        std::vector<int> m_vetoSplitXY;
        Bin2DLeaf* m_leaf;
        unsigned int m_minLeafEntries;
        double m_maxAxisAsymmetry;
        TH2F* m_gridConstraint;

};


#endif
