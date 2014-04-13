

#ifndef BINTREE_H
#define BINTREE_H

#include <iostream>
#include <vector>
#include <algorithm>
#include "TLine.h"
#include "TH2F.h"

class EntryList
{
    public:
        EntryList(int ndim);
        ~EntryList(){};

        void add(const std::vector<double>& values, double weight);

        unsigned int size() const;
        unsigned int effectiveSize() const;
        unsigned int dimension() const;
        double sumOfWeights() const;
        double sumOfWeightsError() const;
        double maxWeight() const;
        double value(unsigned int axis, int entry) const;
        double weight(int entry) const;

        void sort();

        std::pair<EntryList, EntryList> split(unsigned int axis, double cut) const;
        std::pair<int, int> entriesIfSplit(unsigned int axis, double cut) const;

        std::vector<double> percentiles(const std::vector<double>& q, unsigned int axis=0) const;
        double densityGradient(unsigned int axis=0, double q=10.) const;
        void print();


    private:
        unsigned int m_ndim;
        std::vector< std::vector< std::pair<double,int> > > m_sortedValues;
        std::vector< std::vector<int> >  m_sortedPositions;
        std::vector<double> m_weights;
        double m_maxWeight;
        double m_sumOfWeights;
        double m_sumOfWeightsError;
};


class BinLeaf
{
    /* Leaf of a BinTree.
    It encodes a ND bin and stores the entries contained in this bin: (x,y) coordinates + event weight.
    Each leaf is identified with an index.
    */
    public:
        BinLeaf();
        BinLeaf(const std::vector< std::pair<double,double> >& minmax);
        ~BinLeaf();

        bool operator <(const BinLeaf& leaf) const
        {
            return m_index < leaf.index();
        }

        void setBinBoundaries(const std::vector< std::pair<double,double> >& minmax);
        const std::vector< std::pair<double,double> >& getBinBoundaries();
        double getMin(int axis=0);
        double getMax(int axis=0);
        double getWidth(int axis=0);
        double getCenter(int axis=0);
        bool isNeighbor(BinLeaf* leaf);
        unsigned int getNEntries();
        unsigned int effectiveNEntries();
        double getSumOfWeights();
        const EntryList& getEntries();
        //const std::vector< std::vector<double> >& getEntries();
        //const std::vector< double >& getWeights();
        std::vector<double> percentiles(const std::vector<double>& q, unsigned int axis=0);
        double densityGradient(unsigned int axis=0, double q=10.);
        bool inBin(const std::vector<double>& xs);
        bool addEntry(const std::vector<double>& xsi, double wi);
        void setEntries(const EntryList& entries);
        void sortEntries();
        std::vector<TLine*> getBoundaryTLines();

        unsigned int index() const {return m_index;}
        void setIndex(unsigned int index) {m_index = index;}

    private:
        unsigned int m_ndim;
        unsigned int m_index;
        std::vector< std::pair<double,double> > m_binBoundaries;
        EntryList m_entryList;
};


class BinTree
{
    public:
        BinTree(const std::vector< std::pair<double,double> >& minmax, const std::vector< std::vector<double> >& entries, const std::vector< double >& weights);
        ~BinTree();
        void addEntry(const std::vector<double>& xsi, double wi);
        std::vector< std::pair<double,double> > getBinBoundaries();
        double getMin(int axis=0);
        double getMax(int axis=0);
        unsigned int getNEntries();
        double getSumOfWeights();
        //std::vector< std::vector<double> > getEntries();
        //std::vector< double > getWeights();
        BinLeaf* getLeaf(const std::vector<double>& xs);
        std::vector<BinLeaf*> getLeaves();
        std::vector<BinTree*> getTerminalNodes();
        std::vector<BinLeaf*> findNeighborLeaves(BinLeaf* leaf);
        std::vector<BinTree*> getSons() {return m_treeSons;}
        unsigned int getNLeaves();
        double getMinBinWidth(unsigned int axis);
        double getMinEntries();
        double getMaxEntries();
        unsigned int maxLeafIndex();
        unsigned int minLeafEntries() const {return m_minLeafEntries;}
        double maxAxisAsymmetry() const {return m_maxAxisAsymmetry;}
        void build();
        std::vector<TLine*> getBoundaryTLines();
        TH1* fillHistogram();
        std::vector<TH1*> fillWidths(const TH1* widthTemplate=NULL);
        std::vector<TH1*> fillWidthsLowStat(const TH1* widthTemplate=NULL);
        std::vector<TH1*> fillWidthsHighStat(const TH1* widthTemplate=NULL);

        BinLeaf* leaf(){return m_leaf;}
        void setGridConstraint(TH1* gridConstraint);
        void setVetoSplit(unsigned int axis, bool veto){m_vetoSplit[axis] = veto;}
        void setMinLeafEntries(unsigned int minLeafEntries){m_minLeafEntries = minLeafEntries;}
        void setMaxAxisAsymmetry(double maxAxisAsymmetry){m_maxAxisAsymmetry = maxAxisAsymmetry;}
        bool vetoSplit(unsigned int axis){return m_vetoSplit[axis];}

    private:
        std::pair<int,int> entriesIfSplit(double cut, unsigned int axis=0);
        void splitLeaf(double cut, unsigned int maxLeafIndex, unsigned int axis=0);
        void findBestSplit(BinTree*& bestNode, unsigned int& axis, double& gradient);
        void constrainSplit(int axis, double& cut, bool& veto);
        void minimizeLongBins(BinTree* tree, unsigned int axis, double& cut, bool& veto);

        unsigned int m_ndim;
        std::vector<BinTree*> m_treeSons;
        unsigned int m_cutAxis;
        double m_cut;
        std::vector<int> m_vetoSplit;
        BinLeaf* m_leaf;
        unsigned int m_minLeafEntries;
        double m_maxAxisAsymmetry;
        TH1* m_gridConstraint;

};


#endif
