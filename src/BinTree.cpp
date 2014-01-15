#include "BinTree.h"

#include "TAxis.h"
#include "TH3F.h"

#include <iostream>
#include <iomanip>
#include <limits>
#include <map>

using namespace std;


/*****************************************************************/
BinLeaf::BinLeaf()
/*****************************************************************/
{

    m_binBoundaries.push_back(make_pair(0.,1.));
    m_binBoundaries.push_back(make_pair(0.,1.));
    m_index = 0;
    m_ndim = 2;
}

/*****************************************************************/
BinLeaf::~BinLeaf()
/*****************************************************************/
{
    
}

/*****************************************************************/
BinLeaf::BinLeaf(const std::vector< std::pair<double,double> >& minmax)
/*****************************************************************/
{
    unsigned int ndim = minmax.size();
    for(unsigned int axis=0;axis<ndim;axis++)
    {
        if(minmax[axis].first==minmax[axis].second)
        {
            cerr<<"ERROR: BinLeaf.__init__(): Trying to build a bin with zero width ("<<minmax[axis].first<<","<<minmax[axis].second<<")\n";
            exit(0);
        }
    }

    for(unsigned int axis=0;axis<ndim;axis++)
    {
        m_binBoundaries.push_back(make_pair(minmax[axis].first,minmax[axis].second));
    }
    m_index = 0;
    m_ndim = ndim;
}

/*****************************************************************/
void BinLeaf::setBinBoundaries(const std::vector< std::pair<double,double> >& minmax)
/*****************************************************************/
{
    m_binBoundaries.clear();
    for(unsigned int axis=0;axis<m_ndim;axis++)
    {
        m_binBoundaries.push_back(make_pair(minmax[axis].first,minmax[axis].second));
    }
    m_ndim = m_binBoundaries.size();
}

/*****************************************************************/
const std::vector< std::pair<double,double> >& BinLeaf::getBinBoundaries()
/*****************************************************************/
{
    return m_binBoundaries;
}
/*****************************************************************/
double BinLeaf::getMin(int axis)
/*****************************************************************/
{
    return m_binBoundaries[axis].first;
}

/*****************************************************************/
double BinLeaf::getMax(int axis)
/*****************************************************************/
{
    return m_binBoundaries[axis].second;
}


/*****************************************************************/
double BinLeaf::getWidth(int axis)
/*****************************************************************/
{
    return (getMax(axis)-getMin(axis));
}


/*****************************************************************/
double BinLeaf::getCenter(int axis)
/*****************************************************************/
{
    return (getMax(axis)+getMin(axis))/2.;
}



/*****************************************************************/
bool BinLeaf::isNeighbor(BinLeaf* leaf) 
/*****************************************************************/
{
    bool neighbor = false;
    // check if borders are touching
    for(unsigned int axis=0;axis<m_ndim;axis++)
    {
        if( fabs((leaf->getMin(axis)-getMax(axis))/getMax(axis))<1.e-10 )
        {
            for(unsigned int axis2=0;axis2<m_ndim;axis2++)
            {
                if(axis2!=axis)
                {
                    if(leaf->getMax(axis2)>getMin(axis2) && leaf->getMin(axis2)<getMax(axis2) ) neighbor = true;
                }
            }
        }
        if( fabs((leaf->getMax(axis)-getMin(axis))/getMin(axis))<1.e-10 )
        {
            for(unsigned int axis2=0;axis2<m_ndim;axis2++)
            {
                if(axis2!=axis)
                {
                    if(leaf->getMax(axis2)>getMin(axis2) && leaf->getMin(axis2)<getMax(axis2) ) neighbor = true;
                }
            }
        }
    }
    
    return neighbor;
}


/*****************************************************************/
unsigned int BinLeaf::getNEntries()
/*****************************************************************/
{
    return m_entries.size();
}

/*****************************************************************/
double BinLeaf::getSumOfWeights()
/*****************************************************************/
{
    double sumw = 0.;
    for(unsigned int e=0; e<m_weights.size();e++)
    {
        sumw += m_weights[e];
    }
    return sumw;
}

/*****************************************************************/
const std::vector< std::vector<double> >& BinLeaf::getEntries()
/*****************************************************************/
{
    return m_entries;
}

/*****************************************************************/
const std::vector< double >& BinLeaf::getWeights()
/*****************************************************************/
{
    return m_weights;
}

/*****************************************************************/
std::vector<double> BinLeaf::percentiles(const std::vector<double>& q, unsigned int axis)
/*****************************************************************/
{
    if(axis>=m_ndim)
    {
        cerr<<"ERROR: BinLeaf.getPercentiles(): axis>=ndim \n";
        exit(0);
    }
    vector<double> qcopy = q;
    // make sure the quantiles are in increasing order
    sort(qcopy.begin(),qcopy.end());
    // FIXME: don't take into account the entry weights
    vector<double> values;
    for(unsigned int e=0;e<m_entries.size();e++)
    {
        values.push_back(m_entries[e][axis]);
    }
    // sort values
    sort(values.begin(),values.end());
    vector<double> ps;
    //compute quantiles
    double n = 0.;
    double tot = getNEntries();
    double currentq = 0;
    for(unsigned int e=0; e<values.size(); e++)
    {
        n += 1.;
        if(n>=qcopy[currentq]/100.*tot)
        {
            double p = values[e];
            ps.push_back(p);
            currentq++;
            if(currentq>qcopy.size())
            {
                break;
            }
        }
    }
    return ps;
}

/*****************************************************************/
double BinLeaf::densityGradient(unsigned int axis, double q)
/*****************************************************************/
{
    vector<double> qs;
    double qmulti = q;
    while(qmulti<100)
    {
        qs.push_back(qmulti);
        qmulti += q;
    }
    // Filling percentile array
    vector<double> pX = percentiles(qs,axis);
    pX.insert(pX.begin(), m_binBoundaries[axis].first);
    pX.push_back(m_binBoundaries[axis].second);

    //vector<double> localDensities;
    // FIXME: weights are not taken into account here
    unsigned int ntot = getNEntries();
    //double pxtot1 = m_binBoundaries[axis].first;
    //double pxtot2 = m_binBoundaries[axis].second;
    double minDensity = numeric_limits<double>::max();
    double maxDensity = 0.;
    for(unsigned int i=0;i<pX.size()-1;i++)
    {
        double px1 = pX[i];
        double px2 = pX[i+1];
        // Number of entries between two percentiles divided by the distance between them
        double density = ( (float)ntot*(float)q/100. )/(px2-px1);
        if(density<minDensity) minDensity = density;
        if(density>maxDensity) maxDensity = density;
        //localDensities.push_back(density);
    }
    double gradient = fabs(maxDensity-minDensity);
    return gradient;
}

/*****************************************************************/
double BinLeaf::density(double xmin, double xmax, unsigned int axis)
/*****************************************************************/
{
    double nentries = 0.;
    for(unsigned int e=0;e<m_entries.size();e++)
    {
        if(m_entries[e][axis]>=xmin and m_entries[e][axis]<xmax)
        {
            nentries += m_entries[e][2]; // entry weight
        }
    }
    return nentries/(xmax-xmin);
}

/*****************************************************************/
bool BinLeaf::inBin(const std::vector<double>& xs)
/*****************************************************************/
{
    if(xs.size()!=m_ndim)
    {
        cout<<"WARNING: BinLeaf::inBin(): number of dimensions doesn't match\n";
        return false;
    }
    bool isInBin = true;
    for(unsigned int axis=0;axis<m_ndim;axis++)
    {
        if(xs[axis]<getMin(axis) || xs[axis]>getMax(axis)) isInBin = false;
    }
    return isInBin;
}

/*****************************************************************/
bool BinLeaf::addEntry(const std::vector<double>& xsi, double wi)
/*****************************************************************/
{
    if(xsi.size()!=m_ndim)
    {
        cout<<"WARNING: BinLeaf::addEntry(): number of dimensions doesn't match\n";
        return false;
    }
    // Check if it is contained in the bin
    if(inBin(xsi))
    {
        m_entries.push_back(xsi);
        m_weights.push_back(wi);
        return true;
    }
    return false;
}

/*****************************************************************/
std::vector<TLine*> BinLeaf::getBoundaryTLines()
/*****************************************************************/
{
    vector<TLine*> lines;
    if(m_ndim==2)
    {
        TLine* line1 = new TLine(getMin(0), getMin(1), getMin(0), getMax(1));
        TLine* line2 = new TLine(getMin(0), getMax(1), getMax(0), getMax(1));
        TLine* line3 = new TLine(getMax(0), getMax(1), getMax(0), getMin(1));
        TLine* line4 = new TLine(getMax(0), getMin(1), getMin(0), getMin(1));
        lines.push_back(line1);
        lines.push_back(line2);
        lines.push_back(line3);
        lines.push_back(line4);
    }
    return lines;
}


//////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////



/*****************************************************************/
BinTree::BinTree(const std::vector< std::pair<double,double> >& minmax, const std::vector< std::vector<double> >& entries, const std::vector< double >& weights)
/*****************************************************************/
{
    m_treeSons.push_back(NULL);
    m_treeSons.push_back(NULL);
    m_cutAxis = 0;
    m_cut = 0.;
    m_leaf = new BinLeaf(minmax);
    for(unsigned int e=0;e<entries.size();e++)
    {
        m_leaf->addEntry(entries[e], weights[e]);
    }
    m_ndim = minmax.size();
    for(unsigned int axis=0;axis<m_ndim;axis++)
    {
        m_vetoSplit.push_back(false);
    }
    m_minLeafEntries = 200;
    m_maxAxisAsymmetry = 2.;
    m_gridConstraint = NULL;
}

/*****************************************************************/
BinTree::~BinTree()
/*****************************************************************/
{
    //if(m_gridConstraint)
    //{
    //    m_gridConstraint->Delete();
    //    m_gridConstraint = NULL;
    //}
    if(m_leaf)
    {
        delete m_leaf;
    }
    else
    {
        delete m_treeSons[0];
        delete m_treeSons[1];
    }
    
}

void BinTree::setGridConstraint(TH1* gridConstraint) 
{
    //if(m_gridConstraint)
    //{
    //    m_gridConstraint->Delete();
    //    m_gridConstraint = NULL;
    //}
    //m_gridConstraint = (TH2F*)gridConstraint->Clone("gridConstraint");
    //m_gridConstraint->SetDirectory(0);
    m_gridConstraint = gridConstraint;
}

/*****************************************************************/
void BinTree::addEntry(const std::vector<double>& xsi, double wi)
/*****************************************************************/
{
    m_leaf->addEntry(xsi, wi);
}


/*****************************************************************/
std::vector< std::pair<double,double> > BinTree::getBinBoundaries()
/*****************************************************************/
{
    if(m_leaf)
    {
        return m_leaf->getBinBoundaries();
    }
    else
    {
        vector< pair<double, double> > boundaries;
        for(unsigned int axis=0;axis<m_ndim;axis++)
        {
            double min = getMin(axis);
            double max = getMax(axis);
            boundaries.push_back(make_pair(min,max));
        }
        return boundaries;
    }
}

/*****************************************************************/
double BinTree::getMin(int axis)
/*****************************************************************/
{
    if(m_leaf)
    {
        return m_leaf->getMin(axis);
    }
    else
    {
        double xmin1 = m_treeSons[0]->getMin(axis);
        double xmin2 = m_treeSons[1]->getMin(axis);
        double xmin = min(xmin1,xmin2);
        return xmin;
    }
}

/*****************************************************************/
double BinTree::getMax(int axis)
/*****************************************************************/
{
    if(m_leaf)
    {
        return m_leaf->getMax(axis);
    }
    else
    {
        double xmax1 = m_treeSons[0]->getMax(axis);
        double xmax2 = m_treeSons[1]->getMax(axis);
        double xmax = max(xmax1,xmax2);
        return xmax;
    }
}



/*****************************************************************/
unsigned int BinTree::getNEntries()
/*****************************************************************/
{
    if(m_leaf)
    {
        return m_leaf->getNEntries();
    }
    else
    {
        unsigned int nentries = m_treeSons[0]->getNEntries();
        nentries += m_treeSons[1]->getNEntries();
        return nentries;
    }
}


/*****************************************************************/
double BinTree::getSumOfWeights()
/*****************************************************************/
{
    if(m_leaf)
    {
        return m_leaf->getSumOfWeights();
    }
    else
    {
        double sumw = m_treeSons[0]->getSumOfWeights();
        sumw += m_treeSons[1]->getSumOfWeights();
        return sumw;
    }
}


/*****************************************************************/
std::vector< std::vector<double> > BinTree::getEntries()
/*****************************************************************/
{
    if(m_leaf)
    {
        return m_leaf->getEntries();
    }
    else
    {
        vector< vector<double> > entries;
        entries.insert(entries.end(),m_treeSons[0]->getEntries().begin(),m_treeSons[0]->getEntries().end());
        entries.insert(entries.end(),m_treeSons[1]->getEntries().begin(),m_treeSons[1]->getEntries().end());
        return entries;
    }
}

/*****************************************************************/
BinLeaf* BinTree::getLeaf(const std::vector<double>& xs)
/*****************************************************************/
{
    if(m_leaf)
    {
        if(m_leaf->inBin(xs))
        {
            return m_leaf;
        }
        else
        {
            return NULL;
        }
    }
    if(xs[m_cutAxis]<m_cut)
    {
        return m_treeSons[0]->getLeaf(xs);
    }
    else
    {
        return m_treeSons[1]->getLeaf(xs);
    }
}


/*****************************************************************/
std::vector<BinLeaf*> BinTree::getLeaves()
/*****************************************************************/
{
    vector<BinLeaf*> leaves;
    if(m_leaf)
    {
        leaves.push_back(m_leaf);
        return leaves;
    }
    else
    {
        vector<BinLeaf*> leaves1 = m_treeSons[0]->getLeaves();
        vector<BinLeaf*> leaves2 = m_treeSons[1]->getLeaves();
        leaves.insert(leaves.end(),leaves1.begin(), leaves1.end());
        leaves.insert(leaves.end(),leaves2.begin(), leaves2.end());
        return leaves;
    }
}


/*****************************************************************/
std::vector<BinTree*> BinTree::getTerminalNodes()
/*****************************************************************/
{
    vector<BinTree*> nodes;
    if(m_leaf)
    {
        nodes.push_back(this);
        return nodes;
    }
    else
    {
        vector<BinTree*> nodes1 = m_treeSons[0]->getTerminalNodes();
        vector<BinTree*> nodes2 = m_treeSons[1]->getTerminalNodes();
        nodes.insert(nodes.end(),nodes1.begin(), nodes1.end());
        nodes.insert(nodes.end(),nodes2.begin(), nodes2.end());
        return nodes;
    }
}

/*****************************************************************/
std::vector<BinLeaf*> BinTree::findNeighborLeaves(BinLeaf* leaf)
/*****************************************************************/
{
    vector<BinLeaf*> neighborLeaves;
    vector<BinLeaf*> allleaves = getLeaves();
    vector<BinLeaf*>::iterator it = allleaves.begin();
    vector<BinLeaf*>::iterator itE = allleaves.end();
    for(;it!=itE;++it)
    {
        bool isNeighbor = (*it)->isNeighbor(leaf);
        if(isNeighbor)
        {
            neighborLeaves.push_back(*it);
        }
    }
    return neighborLeaves;
}


/*****************************************************************/
unsigned int BinTree::getNLeaves()
/*****************************************************************/
{
    return getLeaves().size();
}

/*****************************************************************/
double BinTree::getMinBinWidth(unsigned int axis)
/*****************************************************************/
{
    if(m_leaf)
    {
        return m_leaf->getWidth(axis);
    }
    else
    {
        double width1 = m_treeSons[0]->getMinBinWidth(axis);
        double width2 = m_treeSons[1]->getMinBinWidth(axis);
        double width = min(width1,width2);
        return width;
    }
}

/*****************************************************************/
double BinTree::getMinEntries()
/*****************************************************************/
{
    if(m_leaf)
    {
        return m_leaf->getNEntries();
    }
    else
    {
        double nentries1 = m_treeSons[0]->getMinEntries();
        double nentries2 = m_treeSons[1]->getMinEntries();
        double nentries = min(nentries1,nentries2);
        return nentries;
    }
}

/*****************************************************************/
double BinTree::getMaxEntries()
/*****************************************************************/
{
    if(m_leaf)
    {
        return m_leaf->getNEntries();
    }
    else
    {
        double nentries1 = m_treeSons[0]->getMaxEntries();
        double nentries2 = m_treeSons[1]->getMaxEntries();
        double nentries = max(nentries1,nentries2);
        return nentries;
    }
}


/*****************************************************************/
unsigned int BinTree::maxLeafIndex()
/*****************************************************************/
{
    if(m_leaf)
    {
        return m_leaf->index();
    }
    else
    {
        return max(m_treeSons[0]->maxLeafIndex(), m_treeSons[1]->maxLeafIndex());
    }
}


/*****************************************************************/
void BinTree::splitLeaf(double cut, unsigned int maxLeafIndex, unsigned int axis)
/*****************************************************************/
{
        if(!m_leaf)
        {
            cerr<<"ERROR: BinTree.split(): This method can only be applied on terminal nodes\n";
            exit(1);
        }
        // Cannot split the bin if the cut value is not contained within the bin boundaries
        if(cut<=getBinBoundaries()[axis].first || cut>=getBinBoundaries()[axis].second)
        {
            cerr<<"ERROR: BinTree.split(): Trying to split a node outside its bin boundaries. "<<cut<<"!=("<<getBinBoundaries()[axis].first<<","<<getBinBoundaries()[axis].second<<")";
            exit(1);
        }
        // Create two neighbour bins
        m_cutAxis = axis;
        m_cut = cut;
        vector< pair<double,double> > boundaries1 = m_leaf->getBinBoundaries();
        vector< pair<double,double> > boundaries2 = m_leaf->getBinBoundaries();
        boundaries1[axis].second = cut;
        boundaries2[axis].first = cut;
        vector< vector<double> > emptyEntries;
        vector< double > emptyWeights;
        m_treeSons[0] = new BinTree(boundaries1, emptyEntries, emptyWeights);
        m_treeSons[1] = new BinTree(boundaries2, emptyEntries, emptyWeights);
        m_treeSons[0]->setMinLeafEntries(m_minLeafEntries);
        m_treeSons[1]->setMinLeafEntries(m_minLeafEntries);
        m_treeSons[0]->setMaxAxisAsymmetry(m_maxAxisAsymmetry);
        m_treeSons[1]->setMaxAxisAsymmetry(m_maxAxisAsymmetry);
        // Fill the two leaves that have just been created with entries of the parent node
        vector< vector<double> > entries = m_leaf->getEntries();
        vector<double> weights = m_leaf->getWeights();
        for(unsigned int e=0;e<entries.size();e++)
        {
            m_treeSons[0]->leaf()->addEntry(entries[e],weights[e]);
            m_treeSons[1]->leaf()->addEntry(entries[e],weights[e]);
        }
        // Assign new leaf indices
        m_treeSons[0]->leaf()->setIndex(maxLeafIndex+1);
        m_treeSons[1]->leaf()->setIndex(maxLeafIndex+2);
        // Set grid constraint
        m_treeSons[0]->setGridConstraint(m_gridConstraint);
        m_treeSons[1]->setGridConstraint(m_gridConstraint);
        // Finally destroy the old leaf. The node is not a terminal node anymore
        delete m_leaf;
        m_leaf = NULL;
}


/*****************************************************************/
void BinTree::findBestSplit(BinTree*& bestNode, unsigned int& axis, double& gradient)
/*****************************************************************/
{
    bestNode = NULL;
    //unsigned int bestAxis = 0;
    // If the node is terminal, only the axis needs to be chosen
    if(m_leaf)
    {
        // Don't split if the bin contains less than 2 times the minimum number of entries
        // FIXME: do we want the stopping condition on Nentries or sum of weights? For the moment it is Nentries
        if(getNEntries()<2.*m_minLeafEntries)
        {
            bestNode = NULL;
            axis = 0;
            gradient = 0.;
            return;
        }
        // Best node is self since this is a terminal node
        bestNode = this;
        // Compute the density gradients along the axis
        // The best axis is the one with the largest gradient
        double maxgrad = 0.;
        int bestAxis = -1;
        for(unsigned int ax=0;ax<m_ndim;ax++)
        {
            double grad = m_leaf->densityGradient(ax);
            if(grad>maxgrad && !m_vetoSplit[ax]) // best gradient and no veto
            {
                maxgrad = grad;
                bestAxis = (int)ax;
            }
        }
        if(bestAxis==-1 || maxgrad==0)
        {
            bestNode = NULL;
            axis = 0;
            gradient = 0.;
            return;
        }
        axis = bestAxis;
        gradient = maxgrad;
        return;
    }
    // If the node is non-terminal, look inside the two sons
    else
    {
        BinTree* bestNode1 = NULL;
        BinTree* bestNode2 = NULL;
        unsigned int bestAxis1 = 0;
        unsigned int bestAxis2 = 0;
        double bestGrad1 = 0.;
        double bestGrad2 = 0.;
        m_treeSons[0]->findBestSplit(bestNode1, bestAxis1, bestGrad1);
        m_treeSons[1]->findBestSplit(bestNode2, bestAxis2, bestGrad2);
        // Look at the largest gradient
        if(bestGrad1>bestGrad2)
        {
            bestNode = bestNode1;
            axis = bestAxis1;
            gradient = bestGrad1;
            return;
        }
        else
        {
            bestNode = bestNode2;
            axis = bestAxis2;
            gradient = bestGrad2;
            return;
        }
    }
}


/*****************************************************************/
void BinTree::constrainSplit(int axis, double& cut, bool& veto)
/*****************************************************************/
{
    if(m_gridConstraint && !m_vetoSplit[axis])
    {
        TAxis* gridAxis = NULL;
        if(axis==0)
        {
            gridAxis = m_gridConstraint->GetXaxis();
        }
        else if(axis==1)
        {
            gridAxis = m_gridConstraint->GetYaxis();
        }
        else if(axis==2)
        {
            gridAxis = m_gridConstraint->GetZaxis();
        }
        else
        {
            cerr<<"ERROR: BinTree::constrainSplit(): Cannot use grid constrain for more than 3D\n";
            exit(1);
        }
        // Find the closest grid constraint for the cut
        // And modify the cut according to this constraint
        int b   = gridAxis->FindBin(cut);
        double low = gridAxis->GetBinLowEdge(b);
        double up  = gridAxis->GetBinUpEdge(b);
        if(fabs(up-cut)<fabs(cut-low))
        {
            cut = up;
            // If the constrained cut is outside the bin boundaries, try the other grid constraint
            if(cut>=getBinBoundaries()[axis].second)
            {
                cut = low;
            }
        }
        else
        {
            cut = low;
            // If the constrained cut is outside the bin boundaries, try the other grid constraint
            if(cut<=getBinBoundaries()[axis].first)
            {
                cut = up;
            }
        }
        //  If the constrained cut is still outside the bin boundaries, veto this bin and axis
        if(cut<=getBinBoundaries()[axis].first || cut>=getBinBoundaries()[axis].second)
        {
            m_vetoSplit[axis] = true;
        }
    }
    veto = m_vetoSplit[axis];
}


/*****************************************************************/
void BinTree::minimizeLongBins(BinTree* tree, unsigned int axis, double& cut, bool& veto)
/*****************************************************************/
{
    if(!tree->vetoSplit(axis))
    {
        vector< pair<double,double> > binBoundaries = tree->getBinBoundaries();
        // this is supposed to be the root tree
        vector< pair<double,double> > fullBoundaries = getBinBoundaries();
        vector<double> fullLengths;
        vector<double> binRelLengths;
        for(unsigned int ax=0;ax<m_ndim;ax++)
        {
            fullLengths.push_back( fullBoundaries[ax].second-fullBoundaries[ax].first );
            binRelLengths.push_back( (binBoundaries[ax].second-binBoundaries[ax].first)/fullLengths[ax] );
        }
        //double binXRelLength = (binBoundaries[axis].second-binBoundaries[axis].first)/fullXLength;
        double cutRelDistance1 = (cut-binBoundaries[axis].first)/fullLengths[axis];
        double cutRelDistance2 = (binBoundaries[axis].second-cut)/fullLengths[axis];
        //unsigned int maxLengthAxis = 0;
        double maxRelLength = 0.;
        for(unsigned int ax=0;ax<m_ndim;ax++)
        {
            if(ax!=axis && binRelLengths[ax]>maxRelLength)
            {
                maxRelLength = binRelLengths[ax];
                //maxLengthAxis = ax;
            }
        }
        //cerr<<" axis="<<axis<<", cut="<<cut<<"\n";
        //cerr<<" cutreldistance = ("<<cutRelDistance1<<","<<cutRelDistance2<<") maxRelLength="<<maxRelLength<<"\n";
        if(cutRelDistance1<cutRelDistance2)
        {
            if(m_maxAxisAsymmetry*cutRelDistance1<maxRelLength)
            {
                cut = maxRelLength/m_maxAxisAsymmetry*fullLengths[axis]+binBoundaries[axis].first;
                cutRelDistance2 = (binBoundaries[axis].second-cut)/fullLengths[axis];
                //cerr<<" 1<max -> cut="<<cut<<"\n";
                if(cut>=binBoundaries[axis].second || m_maxAxisAsymmetry*cutRelDistance2<maxRelLength)
                {
                    tree->setVetoSplit(axis, true);
                    //cerr<<" veto\n";
                }
            }
        }
        else// cutRelDistance2<cutRelDistance1:
        {
            if(m_maxAxisAsymmetry*cutRelDistance2<maxRelLength)
            {
                cut = binBoundaries[axis].second-maxRelLength/m_maxAxisAsymmetry*fullLengths[axis];
                cutRelDistance1 = (cut-binBoundaries[axis].first)/fullLengths[axis];
                //cerr<<" 2<max -> cut="<<cut<<"\n";
                if(cut<=binBoundaries[axis].first || m_maxAxisAsymmetry*cutRelDistance1<maxRelLength)
                {
                    tree->setVetoSplit(axis,true);
                    //cerr<<" veto\n";
                }
            }
        }
    }
    veto = tree->vetoSplit(axis);
}


/*****************************************************************/
void BinTree::build()
/*****************************************************************/
{

    // If the tree already contains too small number of entries, it does nothing
    if(getNEntries()<2.*m_minLeafEntries)
    {
        cout<<"[WARNING] Total number of entries = "<<getNEntries()<<" < 2 x "<<m_minLeafEntries<<". The procedure stops with one single bin\n";
        cout<<"[WARNING]   You'll have to reduce the minimum number of entries per bin if you want to have more than one bin.\n";
        return;
    }
    // Start with first splitting
    BinTree* tree = NULL;
    unsigned int axis = 0;
    double grad = 0.;
    vector<double> perc50;
    perc50.push_back(50.);
    findBestSplit(tree, axis, grad);

    double cut = tree->leaf()->percentiles(perc50,axis)[0];
    // Modify cut according to grid constraints
    bool veto = false;
    tree->constrainSplit(axis, cut, veto);
    if(!veto)
    {
        tree->splitLeaf(cut, maxLeafIndex(), axis);
    }
    int nsplits = 1;
    //int totalEntries = getNEntries();
    //int previousMaxEntries = totalEntries;
    // Split until it is not possible to split (too small number of entries, or vetoed bins)
    while(tree)
    {
        findBestSplit(tree, axis, grad);
        // This is the end
        if(!tree)
        {
            break;
        }
        cut = tree->leaf()->percentiles(perc50,axis)[0];
        veto = false;
        //cerr<<"Bin ("<<tree->getMin(0)<<","<<tree->getMax(0)<<")("<<tree->getMin(1)<<","<<tree->getMax(1)<<")("<<tree->getMin(2)<<","<<tree->getMax(2)<<")\n";
        minimizeLongBins(tree, axis, cut,veto);
        //cerr<<" axis="<<axis<<", cut="<<cut<<"\n";
        if(!veto)
        {
            // Modify cut according to grid constraints
            tree->constrainSplit(axis, cut, veto);
            if(!veto)
            {
                tree->splitLeaf(cut, maxLeafIndex(), axis);
                nsplits += 1;
                cout<<"[INFO]   Number of bins = "<<getNLeaves()<<"\r"<<flush;
            }
        }
        //int maxEntries = getMaxEntries();

        //if( maxEntries!=previousMaxEntries)
        //{
        //    double n = (double)totalEntries/(double)m_minLeafEntries;
        //    double ratio  =  min( 1.-log((double)maxEntries/(double)m_minLeafEntries)/log(n),1.);
        //    int   c      =  ratio * 50;
        //    cout << "[INFO]   "<< setw(3) << (int)(ratio*100) << "% [";
        //    for (int x=0; x<c; x++) cout << "=";
        //    for (int x=c; x<50; x++) cout << " ";
        //    cout << "]\r" << flush;
        //    previousMaxEntries = maxEntries;
        //}

    }


    // Look if some bins are more than 50% empty. If it is the case, the empty part is separated (including one event)
    // from the part of the bin that contains all the entries
    bool finish = false;
    vector<double> perc0;
    vector<double> perc100;
    perc0.push_back(0.);
    perc100.push_back(100.);
    while(!finish)
    {
        vector<BinTree*> terminalNodes = getTerminalNodes();
        for(unsigned int i=0;i<terminalNodes.size();i++)
        {
            BinTree* node = terminalNodes[i];
            BinLeaf* leaf = node->leaf();
            if(leaf->getNEntries()<=1)
            {
                continue;
            }
            vector< pair<double,double> > emptyFractions;
            vector< pair<double,double> > firstLastPoints;
            double maxEmptyFraction = 0.;
            int maxEmptyFractionAxis = 0;
            int maxEmptyFractionSide = 0;
            for(unsigned int axis=0;axis<m_ndim;axis++)
            {
                double firstPoint = leaf->percentiles(perc0,axis)[0];
                double lastPoint  = leaf->percentiles(perc100,axis)[0];
                firstLastPoints.push_back( make_pair(firstPoint, lastPoint) );
                double emptyFractionDown = (firstPoint-leaf->getMin(axis))/(leaf->getMax(axis)-leaf->getMin(axis));
                double emptyFractionUp   = (leaf->getMax(axis)-lastPoint)/(leaf->getMax(axis)-leaf->getMin(axis));
                emptyFractions.push_back( make_pair(emptyFractionDown, emptyFractionUp) );
                if(emptyFractionDown>maxEmptyFraction)
                {
                    maxEmptyFraction = emptyFractionDown;
                    maxEmptyFractionAxis = axis;
                    maxEmptyFractionSide = 0;
                }
                if(emptyFractionUp>maxEmptyFraction)
                {
                    maxEmptyFraction = emptyFractionUp;
                    maxEmptyFractionAxis = axis;
                    maxEmptyFractionSide = 1;
                }
            }
            // If the maximum empty fraction is more than 50%, then the bin is split
            if(maxEmptyFraction>0.5)
            {
                double cut = 0.;
                if(maxEmptyFractionSide==0) cut = firstLastPoints[maxEmptyFractionAxis].first;
                else if(maxEmptyFractionSide==1) cut = firstLastPoints[maxEmptyFractionAxis].second;

                //veto,cut = self.minimizeLongBins(node, cut, axis)
                bool veto = false;
                node->constrainSplit(maxEmptyFractionAxis,cut,veto);
                if(!veto)
                {
                    node->splitLeaf(cut, maxLeafIndex(), maxEmptyFractionAxis);
                }
                nsplits += 1;
            }
            else
            {
                finish = true;
            }
        }
    }
    //cout << "[INFO]   "<< setw(3) << 100 << "% [";
    //for (int x=0; x<50; x++) cout << "=";
    //cout << "]" << endl;
}

/*****************************************************************/
std::vector<TLine*> BinTree::getBoundaryTLines()
/*****************************************************************/
{
    if(m_leaf)
    {
        return m_leaf->getBoundaryTLines();
    }
    else
    {
        vector<TLine*> lines;
        vector<TLine*> lines1 = m_treeSons[0]->getBoundaryTLines();
        vector<TLine*> lines2 = m_treeSons[1]->getBoundaryTLines();
        
        lines.insert(lines.end(),lines1.begin(), lines1.end());
        lines.insert(lines.end(),lines2.begin(), lines2.end());
        return lines;
    }
}


/*****************************************************************/
TH1* BinTree::fillHistogram()
/*****************************************************************/
{
        if(!m_gridConstraint)
        {
            cerr<<"ERROR: BinLeaf::fillHistogram(): Trying to fill histogram, but the binning is unknown. Define first the gridConstraint.\n";
            exit(1);
        }
        if(m_ndim>3)
        {
            cerr<<"ERROR: BinLeaf::fillHistogram(): Cannot fill histograms with more than 3 dimensions\n";
            exit(1);
        }
        TH1* histo = (TH1*)m_gridConstraint->Clone("histoFromTree");
        int nbinsx = histo->GetNbinsX();
        int nbinsy = histo->GetNbinsY();
        int nbinsz = histo->GetNbinsZ();
        map<BinLeaf*, vector< vector<int> > > binsInLeaf;
        // First find the list of TH2 bins for each BinLeaf bin
        for(int bx=1;bx<nbinsx+1;bx++) 
        {
            for(int by=1;by<nbinsy+1;by++)
            {
                if(m_ndim==3)
                {
                    TH3F* h3f = dynamic_cast<TH3F*>(histo);
                    for(int bz=1;bz<nbinsz+1;bz++)
                    {
                        h3f->SetBinContent(bx,by,bz,0);
                        h3f->SetBinError(bx,by,bz,0);
                        double x = h3f->GetXaxis()->GetBinCenter(bx);
                        double y = h3f->GetYaxis()->GetBinCenter(by);
                        double z = h3f->GetZaxis()->GetBinCenter(bz);
                        vector<double> point;
                        point.push_back(x);
                        point.push_back(y);
                        point.push_back(z);
                        BinLeaf* leaf = getLeaf(point);
                        if(binsInLeaf.find(leaf)==binsInLeaf.end())
                        {
                            vector< vector<int> > empty;
                            binsInLeaf[leaf] = empty;
                        }
                        vector<int> bin;
                        bin.push_back(bx);
                        bin.push_back(by);
                        bin.push_back(bz);
                        binsInLeaf[leaf].push_back(bin);
                    }
                }
                else if(m_ndim==2)
                {
                    TH2F* h2f = dynamic_cast<TH2F*>(histo);
                    h2f->SetBinContent(bx,by,0);
                    h2f->SetBinError(bx,by,0);
                    double x = h2f->GetXaxis()->GetBinCenter(bx);
                    double y = h2f->GetYaxis()->GetBinCenter(by);
                    vector<double> point;
                    point.push_back(x);
                    point.push_back(y);
                    BinLeaf* leaf = getLeaf(point);
                    if(binsInLeaf.find(leaf)==binsInLeaf.end())
                    {
                        vector< vector<int> > empty;
                        binsInLeaf[leaf] = empty;
                    }
                    vector<int> bin;
                    bin.push_back(bx);
                    bin.push_back(by);
                    binsInLeaf[leaf].push_back(bin);
                }
            }
        }
        // Then all the TH2 bins are filled according to the entries in the BinLeaf bins
        map<BinLeaf*, vector< vector<int> > >::iterator it = binsInLeaf.begin();
        map<BinLeaf*, vector< vector<int> > >::iterator itE = binsInLeaf.end();
        for(;it!=itE;++it)
        {
            BinLeaf* leaf = it->first;
            vector< vector<int> > bins = it->second;
            vector< vector<double> > entries = leaf->getEntries();
            vector< double > weights = leaf->getWeights();
            int nbins = bins.size();
            for(int b=0;b<nbins;b++)
            {
                if(m_ndim==2)
                {
                    TH2F* h2f = dynamic_cast<TH2F*>(histo);
                    int bx = bins[b][0];
                    int by = bins[b][1];
                    double x = h2f->GetXaxis()->GetBinCenter(bx);
                    double y = h2f->GetYaxis()->GetBinCenter(by);
                    for(unsigned int e=0;e<entries.size();e++)
                    {
                        double value = weights[e]/(double)nbins;
                        h2f->Fill(x,y,value);
                    }
                }
                else if(m_ndim==3)
                {
                    TH3F* h3f = dynamic_cast<TH3F*>(histo);
                    int bx = bins[b][0];
                    int by = bins[b][1];
                    int bz = bins[b][2];
                    double x = h3f->GetXaxis()->GetBinCenter(bx);
                    double y = h3f->GetYaxis()->GetBinCenter(by);
                    double z = h3f->GetZaxis()->GetBinCenter(bz);
                    for(unsigned int e=0;e<entries.size();e++)
                    {
                        double value = weights[e]/(double)nbins;
                        h3f->Fill(x,y,z,value);
                    }
                }
            }
        }
        return histo;
}

/*****************************************************************/
vector<TH1*> BinTree::fillWidths(const TH1* widthTemplate)
/*****************************************************************/
{
        if(!widthTemplate && !m_gridConstraint)
        {
            cerr<<"ERROR: BinLeaf::fillWidths(): Trying to fill width, but the binning is unknown. Please give a template histogram of define a gridConstraint.\n";
            exit(1);
        }
        if(m_ndim>3)
        {
            cerr<<"ERROR: BinLeaf::fillWidths(): Cannot fill histograms with more than 3 dimensions\n";
            exit(1);
        }
        vector<TH1*> widths;
        if(m_ndim==2)
        {
            const TH1* gridRef = (widthTemplate ? widthTemplate : m_gridConstraint);
            TH1* hWidthX = (TH1*)gridRef->Clone("widthXFromTree");
            TH1* hWidthY = (TH1*)gridRef->Clone("widthYFromTree");
            int nbinsx = hWidthX->GetNbinsX();
            int nbinsy = hWidthX->GetNbinsY();
            int counter = 0;
            int total = nbinsx*nbinsy;
            for(int bx=1;bx<nbinsx+1;bx++) 
            {
                for(int by=1;by<nbinsy+1;by++)
                {
                    counter++;
                    if(counter % (total/100) == 0)
                    {
                        double ratio  =  (double)counter/(double)total;
                        int   c      =  ratio * 50;
                        cout << "[INFO]   "<< setw(3) << (int)(ratio*100) << "% [";
                        for (int x=0; x<c; x++) cout << "=";
                        for (int x=c; x<50; x++) cout << " ";
                        cout << "]\r" << flush;
                    }
                    double x = hWidthX->GetXaxis()->GetBinCenter(bx);
                    double y = hWidthX->GetYaxis()->GetBinCenter(by);
                    vector<double> point;
                    point.push_back(x);
                    point.push_back(y);
                    //cout<<"Computing width for ("<<x<<","<<y<<")\n";
                    BinLeaf* leaf = getLeaf(point);
                    vector<BinLeaf*> neighborLeaves = findNeighborLeaves(leaf);
                    neighborLeaves.push_back(leaf);
                    vector<BinLeaf*>::iterator itLeaf = neighborLeaves.begin();
                    vector<BinLeaf*>::iterator itELeaf = neighborLeaves.end();
                    double sumw = 0.;
                    double sumwx = 0.;
                    double sumwy = 0.;
                    //cout<<" "<<neighborLeaves.size()<<" neighbor leaves\n";
                    for(;itLeaf!=itELeaf;++itLeaf)
                    {
                        double xi = (*itLeaf)->getCenter(0);
                        double yi = (*itLeaf)->getCenter(1);
                        double dx = fabs(xi-x);
                        double dy = fabs(yi-y);
                        double wxi = (*itLeaf)->getWidth(0);
                        double wyi = (*itLeaf)->getWidth(1);
                        if(dx<0.05*wxi) dx = 0.05*wxi;
                        if(dy<0.05*wyi) dy = 0.05*wyi;
                        double dr2 = dx*dx+dy*dy;
                        double dr = sqrt(dr2);
                        sumw += 1./dr;
                        sumwx += wxi/dr;
                        sumwy += wyi/dr;
                        //cout<<"  Neighbor leaf dx="<<dx<<",dy="<<dy<<",w="<<1./dr2<<"\n";
                    }
                    double widthx = sumwx/sumw;
                    double widthy = sumwy/sumw;
                    //cout<<" Width x="<<widthx<<",y="<<widthy<<"\n";
                    hWidthX->SetBinContent(bx,by,widthx);
                    hWidthY->SetBinContent(bx,by,widthy);
                    hWidthX->SetBinError(bx,by,0.);
                    hWidthY->SetBinError(bx,by,0.);
                }
            }
            widths.push_back(hWidthX);
            widths.push_back(hWidthY);
        }
        else if(m_ndim==3)
        {
            const TH1* gridRef = (widthTemplate ? widthTemplate : m_gridConstraint);
            TH1* hWidthX = (TH1*)gridRef->Clone("widthXFromTree");
            TH1* hWidthY = (TH1*)gridRef->Clone("widthYFromTree");
            TH1* hWidthZ = (TH1*)gridRef->Clone("widthZFromTree");
            int nbinsx = hWidthX->GetNbinsX();
            int nbinsy = hWidthX->GetNbinsY();
            int nbinsz = hWidthX->GetNbinsZ();
            int counter = 0;
            int total = nbinsx*nbinsy*nbinsz;
            for(int bx=1;bx<nbinsx+1;bx++) 
            {
                for(int by=1;by<nbinsy+1;by++)
                {
                    for(int bz=1;bz<nbinsz+1;bz++)
                    {
                        counter++;
                        if(counter % (total/100) == 0)
                        {
                            double ratio  =  (double)counter/(double)total;
                            int   c      =  ratio * 50;
                            cout << "[INFO]   "<< setw(3) << (int)(ratio*100) << "% [";
                            for (int x=0; x<c; x++) cout << "=";
                            for (int x=c; x<50; x++) cout << " ";
                            cout << "]\r" << flush;
                        }
                        double x = hWidthX->GetXaxis()->GetBinCenter(bx);
                        double y = hWidthX->GetYaxis()->GetBinCenter(by);
                        double z = hWidthX->GetZaxis()->GetBinCenter(bz);
                        vector<double> point;
                        point.push_back(x);
                        point.push_back(y);
                        point.push_back(z);
                        //cout<<"Computing width for ("<<x<<","<<y<<")\n";
                        BinLeaf* leaf = getLeaf(point);
                        vector<BinLeaf*> neighborLeaves = findNeighborLeaves(leaf);
                        neighborLeaves.push_back(leaf);
                        vector<BinLeaf*>::iterator itLeaf = neighborLeaves.begin();
                        vector<BinLeaf*>::iterator itELeaf = neighborLeaves.end();
                        double sumw = 0.;
                        double sumwx = 0.;
                        double sumwy = 0.;
                        double sumwz = 0.;
                        //cout<<" "<<neighborLeaves.size()<<" neighbor leaves\n";
                        for(;itLeaf!=itELeaf;++itLeaf)
                        {
                            double xi = (*itLeaf)->getCenter(0);
                            double yi = (*itLeaf)->getCenter(1);
                            double zi = (*itLeaf)->getCenter(2);
                            double dx = fabs(xi-x);
                            double dy = fabs(yi-y);
                            double dz = fabs(zi-z);
                            double wxi = (*itLeaf)->getWidth(0);
                            double wyi = (*itLeaf)->getWidth(1);
                            double wzi = (*itLeaf)->getWidth(2);
                            if(dx<0.05*wxi) dx = 0.05*wxi;
                            if(dy<0.05*wyi) dy = 0.05*wyi;
                            if(dz<0.05*wzi) dz = 0.05*wzi;
                            double dr2 = dx*dx+dy*dy+dz*dz;
                            double dr = sqrt(dr2);
                            sumw += 1./dr;
                            sumwx += wxi/dr;
                            sumwy += wyi/dr;
                            sumwz += wzi/dr;
                            //cout<<"  Neighbor leaf dx="<<dx<<",dy="<<dy<<",w="<<1./dr2<<"\n";
                        }
                        double widthx = sumwx/sumw;
                        double widthy = sumwy/sumw;
                        double widthz = sumwz/sumw;
                        hWidthX->SetBinContent(bx,by,bz,widthx);
                        hWidthY->SetBinContent(bx,by,bz,widthy);
                        hWidthZ->SetBinContent(bx,by,bz,widthz);
                        hWidthX->SetBinError(bx,by,bz,0.);
                        hWidthY->SetBinError(bx,by,bz,0.);
                        hWidthZ->SetBinError(bx,by,bz,0.);
                    }
                }
            }
            widths.push_back(hWidthX);
            widths.push_back(hWidthY);
            widths.push_back(hWidthZ);
        }
        cout << "[INFO]   "<< setw(3) << 100 << "% [";
        for (int x=0; x<50; x++) cout << "=";
        cout << "]" << endl;
        return widths;
}





