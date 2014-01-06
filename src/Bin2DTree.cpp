#include "Bin2DTree.h"

#include "TAxis.h"

#include <iostream>
#include <limits>
#include <map>

using namespace std;

/*****************************************************************/
Bin2DLeaf::Bin2DLeaf()
/*****************************************************************/
{

    m_binBoundaries.push_back(make_pair(0.,1.));
    m_binBoundaries.push_back(make_pair(0.,1.));
    m_index = 0;
}

/*****************************************************************/
Bin2DLeaf::~Bin2DLeaf()
/*****************************************************************/
{
    
}

/*****************************************************************/
Bin2DLeaf::Bin2DLeaf(double xmin, double xmax, double ymin, double ymax)
/*****************************************************************/
{
    if(xmin==xmax || ymin==ymax)
    {
        cerr<<"ERROR: Bin2DLeaf.__init__(): Trying to build a bin with zero width ("<<xmin<<","<<xmax<<"),("<<ymin<<","<<ymax<<")\n";
        exit(0);
    }

    m_binBoundaries.push_back(make_pair(xmin,xmax));
    m_binBoundaries.push_back(make_pair(ymin,ymax));
    m_index = 0;
}

/*****************************************************************/
void Bin2DLeaf::setBinBoundaries(double xmin, double xmax, double ymin, double ymax)
/*****************************************************************/
{
    m_binBoundaries.clear();
    m_binBoundaries.push_back(make_pair(xmin,xmax));
    m_binBoundaries.push_back(make_pair(ymin,ymax));
}

/*****************************************************************/
const std::vector< std::pair<double,double> >& Bin2DLeaf::getBinBoundaries()
/*****************************************************************/
{
    return m_binBoundaries;
}
/*****************************************************************/
double Bin2DLeaf::getXMin()
/*****************************************************************/
{
    return m_binBoundaries[0].first;
}

/*****************************************************************/
double Bin2DLeaf::getXMax()
/*****************************************************************/
{
    return m_binBoundaries[0].second;
}

/*****************************************************************/
double Bin2DLeaf::getYMin()
/*****************************************************************/
{
    return m_binBoundaries[1].first;
}

/*****************************************************************/
double Bin2DLeaf::getYMax()
/*****************************************************************/
{
    return m_binBoundaries[1].second;
}

/*****************************************************************/
double Bin2DLeaf::getXWidth()
/*****************************************************************/
{
    return (getXMax()-getXMin());
}

/*****************************************************************/
double Bin2DLeaf::getYWidth()
/*****************************************************************/
{
    return (getYMax()-getYMin());
}

/*****************************************************************/
double Bin2DLeaf::getXCenter()
/*****************************************************************/
{
    return (getXMax()+getXMin())/2.;
}

/*****************************************************************/
double Bin2DLeaf::getYCenter()
/*****************************************************************/
{
    return (getYMax()+getYMin())/2.;
}


/*****************************************************************/
bool Bin2DLeaf::isNeighbor(Bin2DLeaf* leaf) 
/*****************************************************************/
{
    bool neighbor = false;
    if( fabs((leaf->getXMin()-getXMax())/getXMax())<1.e-10 ) neighbor = true;
    if( fabs((leaf->getXMax()-getXMin())/getXMin())<1.e-10 ) neighbor = true;
    if( fabs((leaf->getYMin()-getYMax())/getYMax())<1.e-10 ) neighbor = true;
    if( fabs((leaf->getYMax()-getYMin())/getYMin())<1.e-10 ) neighbor = true;
    
    return neighbor;
}


/*****************************************************************/
unsigned int Bin2DLeaf::getNEntries()
/*****************************************************************/
{
    return m_entries.size();
}

/*****************************************************************/
double Bin2DLeaf::getSumOfWeights()
/*****************************************************************/
{
    double sumw = 0.;
    for(unsigned int e=0; e<m_entries.size();e++)
    {
        sumw += m_entries[e][2];
    }
    return sumw;
}

/*****************************************************************/
const std::vector< std::vector<double> >& Bin2DLeaf::getEntries()
/*****************************************************************/
{
    return m_entries;
}

/*****************************************************************/
std::vector<double> Bin2DLeaf::percentiles(const std::vector<double>& q, unsigned int axis)
/*****************************************************************/
{
    if(axis>1)
    {
        cerr<<"ERROR: Bin2DLeaf.getPercentiles(): axis>1 \n";
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
double Bin2DLeaf::densityGradient(unsigned int axis, double q)
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
    double pytot1 = m_binBoundaries[(axis+1)%2].first;
    double pytot2 = m_binBoundaries[(axis+1)%2].second;
    double minDensity = numeric_limits<double>::max();
    double maxDensity = 0.;
    for(unsigned int i=0;i<pX.size()-1;i++)
    {
        double px1 = pX[i];
        double px2 = pX[i+1];
        // Number of entries between two percentiles divided by the distance between them
        double density = ( (float)ntot*(float)q/100. )/( (px2-px1)*(pytot2-pytot1) );
        if(density<minDensity) minDensity = density;
        if(density>maxDensity) maxDensity = density;
        //localDensities.push_back(density);
    }
    double gradient = fabs(maxDensity-minDensity);
    return gradient;
}

/*****************************************************************/
double Bin2DLeaf::density(double xmin, double xmax, unsigned int axis)
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
bool Bin2DLeaf::inBin(double xi, double yi)
/*****************************************************************/
{
    return (xi>=getXMin() && xi<getXMax() && yi>=getYMin() && yi<getYMax());
}

/*****************************************************************/
bool Bin2DLeaf::addEntry(double xi, double yi, double wi)
/*****************************************************************/
{
    // First check if it is contained in the bin
    if(inBin(xi,yi))
    {
        vector<double> entry;
        entry.push_back(xi);
        entry.push_back(yi);
        entry.push_back(wi);
        m_entries.push_back(entry);
        return true;
    }
    return false;
}

/*****************************************************************/
std::vector<TLine*> Bin2DLeaf::getBoundaryTLines()
/*****************************************************************/
{
    TLine* line1 = new TLine(getXMin(), getYMin(), getXMin(), getYMax());
    TLine* line2 = new TLine(getXMin(), getYMax(), getXMax(), getYMax());
    TLine* line3 = new TLine(getXMax(), getYMax(), getXMax(), getYMin());
    TLine* line4 = new TLine(getXMax(), getYMin(), getXMin(), getYMin());
    vector<TLine*> lines;
    lines.push_back(line1);
    lines.push_back(line2);
    lines.push_back(line3);
    lines.push_back(line4);
    return lines;
}


//////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////



/*****************************************************************/
Bin2DTree::Bin2DTree(double xmin, double xmax, double ymin, double ymax, const std::vector< std::vector<double> >& entries)
/*****************************************************************/
{
    m_treeSons.push_back(NULL);
    m_treeSons.push_back(NULL);
    m_cutAxis = 0;
    m_cut = 0.;
    m_leaf = new Bin2DLeaf(xmin,xmax,ymin,ymax);
    for(unsigned int e=0;e<entries.size();e++)
    {
        m_leaf->addEntry(entries[e][0],entries[e][1],entries[e][2]);
    }
    m_vetoSplitXY.push_back(false);
    m_vetoSplitXY.push_back(false);
    m_minLeafEntries = 200;
    m_maxAxisAsymmetry = 2.;
    m_gridConstraint = NULL;
}

/*****************************************************************/
Bin2DTree::~Bin2DTree()
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

void Bin2DTree::setGridConstraint(TH2F* gridConstraint) 
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
void Bin2DTree::addEntry(double xi, double yi, double wi)
/*****************************************************************/
{
    m_leaf->addEntry(xi,yi,wi);
}


/*****************************************************************/
std::vector< std::pair<double,double> > Bin2DTree::getBinBoundaries()
/*****************************************************************/
{
    if(m_leaf)
    {
        return m_leaf->getBinBoundaries();
    }
    else
    {
        double xmin = getXMin();
        double xmax = getXMax();
        double ymin = getYMin();
        double ymax = getYMax();
        vector< pair<double, double> > boundaries;
        boundaries.push_back(make_pair(xmin,xmax));
        boundaries.push_back(make_pair(ymin,ymax));
        return boundaries;
    }
}

/*****************************************************************/
double Bin2DTree::getXMin()
/*****************************************************************/
{
    if(m_leaf)
    {
        return m_leaf->getXMin();
    }
    else
    {
        double xmin1 = m_treeSons[0]->getXMin();
        double xmin2 = m_treeSons[1]->getXMin();
        double xmin = min(xmin1,xmin2);
        return xmin;
    }
}

/*****************************************************************/
double Bin2DTree::getXMax()
/*****************************************************************/
{
    if(m_leaf)
    {
        return m_leaf->getXMax();
    }
    else
    {
        double xmax1 = m_treeSons[0]->getXMax();
        double xmax2 = m_treeSons[1]->getXMax();
        double xmax = max(xmax1,xmax2);
        return xmax;
    }
}


/*****************************************************************/
double Bin2DTree::getYMin()
/*****************************************************************/
{
    if(m_leaf)
    {
        return m_leaf->getYMin();
    }
    else
    {
        double ymin1 = m_treeSons[0]->getYMin();
        double ymin2 = m_treeSons[1]->getYMin();
        double ymin = min(ymin1,ymin2);
        return ymin;
    }
}

/*****************************************************************/
double Bin2DTree::getYMax()
/*****************************************************************/
{
    if(m_leaf)
    {
        return m_leaf->getYMax();
    }
    else
    {
        double ymax1 = m_treeSons[0]->getYMax();
        double ymax2 = m_treeSons[1]->getYMax();
        double ymax = max(ymax1,ymax2);
        return ymax;
    }
}


/*****************************************************************/
unsigned int Bin2DTree::getNEntries()
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
double Bin2DTree::getSumOfWeights()
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
std::vector< std::vector<double> > Bin2DTree::getEntries()
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
Bin2DLeaf* Bin2DTree::getLeaf(double xi, double yi)
/*****************************************************************/
{
    vector<double> x;
    x.push_back(xi);
    x.push_back(yi);
    if(m_leaf)
    {
        if(m_leaf->inBin(xi,yi))
        {
            return m_leaf;
        }
        else
        {
            return NULL;
        }
    }
    if(x[m_cutAxis]<m_cut)
    {
        return m_treeSons[0]->getLeaf(xi,yi);
    }
    else
    {
        return m_treeSons[1]->getLeaf(xi,yi);
    }
}


/*****************************************************************/
std::vector<Bin2DLeaf*> Bin2DTree::getLeaves()
/*****************************************************************/
{
    vector<Bin2DLeaf*> leaves;
    if(m_leaf)
    {
        leaves.push_back(m_leaf);
        return leaves;
    }
    else
    {
        vector<Bin2DLeaf*> leaves1 = m_treeSons[0]->getLeaves();
        vector<Bin2DLeaf*> leaves2 = m_treeSons[1]->getLeaves();
        leaves.insert(leaves.end(),leaves1.begin(), leaves1.end());
        leaves.insert(leaves.end(),leaves2.begin(), leaves2.end());
        return leaves;
    }
}


/*****************************************************************/
std::vector<Bin2DTree*> Bin2DTree::getTerminalNodes()
/*****************************************************************/
{
    vector<Bin2DTree*> nodes;
    if(m_leaf)
    {
        nodes.push_back(this);
        return nodes;
    }
    else
    {
        vector<Bin2DTree*> nodes1 = m_treeSons[0]->getTerminalNodes();
        vector<Bin2DTree*> nodes2 = m_treeSons[1]->getTerminalNodes();
        nodes.insert(nodes.end(),nodes1.begin(), nodes1.end());
        nodes.insert(nodes.end(),nodes2.begin(), nodes2.end());
        return nodes;
    }
}

/*****************************************************************/
std::vector<Bin2DLeaf*> Bin2DTree::findNeighborLeaves(Bin2DLeaf* leaf)
/*****************************************************************/
{
    vector<Bin2DLeaf*> neighborLeaves;
    vector<Bin2DLeaf*> allleaves = getLeaves();
    vector<Bin2DLeaf*>::iterator it = allleaves.begin();
    vector<Bin2DLeaf*>::iterator itE = allleaves.end();
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
unsigned int Bin2DTree::getNLeaves()
/*****************************************************************/
{
    return getLeaves().size();
}


/*****************************************************************/
unsigned int Bin2DTree::maxLeafIndex()
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
void Bin2DTree::splitLeaf(double cut, unsigned int maxLeafIndex, unsigned int axis)
/*****************************************************************/
{
        if(!m_leaf)
        {
            cerr<<"ERROR: Bin2DTree.split(): This method can only be applied on terminal nodes\n";
            exit(1);
        }
        // Cannot split the bin if the cut value is not contained within the bin boundaries
        if(cut<=getBinBoundaries()[axis].first || cut>=getBinBoundaries()[axis].second)
        {
            cerr<<"ERROR: Bin2DTree.split(): Trying to split a node outside its bin boundaries. "<<cut<<"!=("<<getBinBoundaries()[axis].first<<","<<getBinBoundaries()[axis].second<<")";
        }
        // Create two neighbouring bins
        m_cutAxis = axis;
        m_cut = cut;
        vector< pair<double,double> > boundaries1 = m_leaf->getBinBoundaries();
        vector< pair<double,double> > boundaries2 = m_leaf->getBinBoundaries();
        boundaries1[axis].second = cut;
        boundaries2[axis].first = cut;
        vector< vector<double> > emptyEntries;
        m_treeSons[0] = new Bin2DTree(boundaries1[0].first,boundaries1[0].second,boundaries1[1].first,boundaries1[1].second, emptyEntries);
        m_treeSons[1] = new Bin2DTree(boundaries2[0].first,boundaries2[0].second,boundaries2[1].first,boundaries2[1].second, emptyEntries);
        m_treeSons[0]->setMinLeafEntries(m_minLeafEntries);
        m_treeSons[1]->setMinLeafEntries(m_minLeafEntries);
        m_treeSons[0]->setMaxAxisAsymmetry(m_maxAxisAsymmetry);
        m_treeSons[1]->setMaxAxisAsymmetry(m_maxAxisAsymmetry);
        // Fill the two leaves that have just been created with entries of the parent node
        vector< vector<double> > entries = m_leaf->getEntries();
        for(unsigned int e=0;e<entries.size();e++)
        {
            m_treeSons[0]->leaf()->addEntry(entries[e][0],entries[e][1],entries[e][2]);
            m_treeSons[1]->leaf()->addEntry(entries[e][0],entries[e][1],entries[e][2]);
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
void Bin2DTree::findBestSplit(Bin2DTree*& bestNode, unsigned int& axis, double& gradient)
/*****************************************************************/
{
    bestNode = NULL;
    unsigned int bestAxis = 0;
    // If the node is terminal, only the axis needs to be chosen
    if(m_leaf)
    {
        // Don't split if the bin contains less than 2 times the minimum number of entries
        // FIXME: do we want the stopping condition on Nentries or sum of weights?
        if(getNEntries()<2.*m_minLeafEntries)
        {
            bestNode = NULL;
            axis = 0;
            gradient = 0.;
            return;
        }
        // Best node is self since this is a terminal node
        bestNode = this;
        // Compute the density gradients along the two axis
        vector<double> grads;
        grads.push_back(m_leaf->densityGradient(0));
        grads.push_back(m_leaf->densityGradient(1));
        // The best axis is the one with the largest gradient
        if(grads[0]>=grads[1])
        {
            bestAxis = 0;
        }
        else
        {
            bestAxis = 1;
        }
        // If there is a veto on this node and on the chosen axis, the other axis is the only possible choice
        if(m_vetoSplitXY[bestAxis])
        {
            bestAxis = (bestAxis+1)%2;
            // If there is a veto on both axis, then this bin cannot be split
            if(m_vetoSplitXY[bestAxis])
            {
                bestNode = NULL;
                axis = 0;
                gradient = 0.;
                return;
            }
        }
        if(grads[bestAxis]==0)
        {
            bestNode = NULL;
            axis = 0;
            gradient = 0.;
            return;
        }

        axis = bestAxis;
        gradient = grads[bestAxis];
        return;
    }
    // If the node is non-terminal, look inside the two sons
    else
    {
        Bin2DTree* bestNode1 = NULL;
        Bin2DTree* bestNode2 = NULL;
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
void Bin2DTree::constrainSplit(int axis, double& cut, bool& veto)
/*****************************************************************/
{
    if(m_gridConstraint && !m_vetoSplitXY[axis])
    {
        TAxis* gridAxis = NULL;
        if(axis==0)
        {
            gridAxis = m_gridConstraint->GetXaxis();
        }
        else
        {
            gridAxis = m_gridConstraint->GetYaxis();
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
            m_vetoSplitXY[axis] = true;
        }
    }
    veto = m_vetoSplitXY[axis];
}


/*****************************************************************/
void Bin2DTree::minimizeLongBins(Bin2DTree* tree, unsigned int axis, double& cut, bool& veto)
/*****************************************************************/
{
    if(!m_vetoSplitXY[axis])
    {
        vector< pair<double,double> > binBoundaries = tree->getBinBoundaries();
        // this is supposed to be the root tree
        vector< pair<double,double> > fullBoundaries = getBinBoundaries();
        double fullXLength = fullBoundaries[axis].second-fullBoundaries[axis].first;
        double fullYLength = fullBoundaries[(axis+1)%2].second-fullBoundaries[(axis+1)%2].first;
        //double binXRelLength = (binBoundaries[axis].second-binBoundaries[axis].first)/fullXLength;
        double binYRelLength = (binBoundaries[(axis+1)%2].second-binBoundaries[(axis+1)%2].first)/fullYLength;
        double cutRelDistance1 = (cut-binBoundaries[axis].first)/fullXLength;
        double cutRelDistance2 = (binBoundaries[axis].second-cut)/fullXLength;
        if(cutRelDistance1<cutRelDistance2)
        {
            if(m_maxAxisAsymmetry*cutRelDistance1<binYRelLength)
            {
                cut = binYRelLength/m_maxAxisAsymmetry*fullXLength+binBoundaries[axis].first;
                cutRelDistance2 = (binBoundaries[axis].second-cut)/fullXLength;
                if(cut>=binBoundaries[axis].second || m_maxAxisAsymmetry*cutRelDistance2<binYRelLength)
                {
                    tree->setVetoSplitXY(axis, true);
                }
            }
        }
        else// cutRelDistance2<cutRelDistance1:
        {
            if(m_maxAxisAsymmetry*cutRelDistance2<binYRelLength)
            {
                cut = binBoundaries[axis].second-binYRelLength/m_maxAxisAsymmetry*fullXLength;
                cutRelDistance1 = (cut-binBoundaries[axis].first)/fullXLength;
                if(cut<=binBoundaries[axis].first || m_maxAxisAsymmetry*cutRelDistance1<binYRelLength)
                {
                    tree->setVetoSplitXY(axis,true);
                }
            }
        }
    }
    veto = tree->vetoSplitXY(axis);
}


/*****************************************************************/
void Bin2DTree::build()
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
    Bin2DTree* tree = NULL;
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
        minimizeLongBins(tree, axis, cut,veto);
        // Modify cut according to grid constraints
        tree->constrainSplit(axis, cut, veto);
        if(!veto)
        {
            tree->splitLeaf(cut, maxLeafIndex(), axis);
        }
        nsplits += 1;
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
        vector<Bin2DTree*> terminalNodes = getTerminalNodes();
        for(unsigned int i=0;i<terminalNodes.size();i++)
        {
            Bin2DTree* node = terminalNodes[i];
            Bin2DLeaf* leaf = node->leaf();
            if(leaf->getNEntries()<=1)
            {
                continue;
            }
            // look for the x and y values of the first and last entries
            double firstPointX = leaf->percentiles(perc0,0)[0];
            double firstPointY = leaf->percentiles(perc0,1)[0];
            double lastPointX = leaf->percentiles(perc100,0)[0];
            double lastPointY = leaf->percentiles(perc100,1)[0];
            // compute the empty fractions along the two axis
            double emptyFractionLeft = (firstPointX-leaf->getXMin())/(leaf->getXMax()-leaf->getXMin());
            double emptyFractionRight = (leaf->getXMax()-lastPointX)/(leaf->getXMax()-leaf->getXMin());
            double emptyFractionBottom = (firstPointY-leaf->getYMin())/(leaf->getYMax()-leaf->getYMin());
            double emptyFractionTop = (leaf->getYMax()-lastPointY)/(leaf->getYMax()-leaf->getYMin());
            // Take the maximum empty fraction
            double maxEmptyFraction = max(emptyFractionLeft, max(emptyFractionRight, max(emptyFractionBottom,emptyFractionTop)));
            // If the maximum empty fraction is more tha 50%, then the bin is split
            if(maxEmptyFraction>0.5)
            {
                double cut = 0;
                unsigned int axis = 0;
                if(maxEmptyFraction==emptyFractionLeft)
                {
                    cut = firstPointX;
                    axis = 0;
                }
                else if( maxEmptyFraction==emptyFractionRight)
                {
                    cut = lastPointX;
                    axis = 0;
                }
                else if( maxEmptyFraction==emptyFractionBottom)
                {
                    cut = firstPointY;
                    axis = 1;
                }
                else //maxEmptyFraction==emptyFractionTop:
                {
                    cut = lastPointY;
                    axis = 1;
                }
                //veto,cut = self.minimizeLongBins(node, cut, axis)
                bool veto = false;
                node->constrainSplit(axis,cut,veto);
                if(!veto)
                {
                    node->splitLeaf(cut, maxLeafIndex(), axis);
                }
                nsplits += 1;
            }
            else
            {
                finish = true;
            }
        }
    }
}

/*****************************************************************/
std::vector<TLine*> Bin2DTree::getBoundaryTLines()
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
TH2F* Bin2DTree::fillHistogram()
/*****************************************************************/
{
        if(!m_gridConstraint)
        {
            cerr<<"ERROR: Bin2DLeaf.fillHistogram(): Trying to fill histogram, but the binning is unknown. Define first the gridConstraint.\n";
            exit(1);
        }
        TH2F* histo = (TH2F*)m_gridConstraint->Clone("histoFromTree");
        int nbinsx = histo->GetNbinsX();
        int nbinsy = histo->GetNbinsY();
        map<Bin2DLeaf*, vector< pair<int,int> > > binsInLeaf;
        // First find the list of TH2 bins for each Bin2DLeaf bin
        for(int bx=1;bx<nbinsx+1;bx++) 
        {
            for(int by=1;by<nbinsy+1;by++)
            {
                histo->SetBinContent(bx,by,0);
                histo->SetBinError(bx,by,0);
                double x = histo->GetXaxis()->GetBinCenter(bx);
                double y = histo->GetYaxis()->GetBinCenter(by);
                Bin2DLeaf* leaf = getLeaf(x,y);
                if(binsInLeaf.find(leaf)==binsInLeaf.end())
                {
                    vector< pair<int,int> > empty;
                    binsInLeaf[leaf] = empty;
                }
                binsInLeaf[leaf].push_back(make_pair(bx,by));
            }
        }
        // Then all the TH2 bins are filled according to the entries in the Bin2DLeaf bins
        map<Bin2DLeaf*, vector< pair<int,int> > >::iterator it = binsInLeaf.begin();
        map<Bin2DLeaf*, vector< pair<int,int> > >::iterator itE = binsInLeaf.end();
        for(;it!=itE;++it)
        {
            Bin2DLeaf* leaf = it->first;
            vector< pair<int,int> > bins = it->second;
            vector< vector<double> > entries = leaf->getEntries();
            int nbins = bins.size();
            for(int b=0;b<nbins;b++)
            {
                int bx = bins[b].first;
                int by = bins[b].second;
                double x = histo->GetXaxis()->GetBinCenter(bx);
                double y = histo->GetYaxis()->GetBinCenter(by);
                for(unsigned int e=0;e<entries.size();e++)
                {
                    double value = entries[e][2]/(double)nbins;
                    histo->Fill(x,y,value);
                }
            }
        }
        return histo;
}

/*****************************************************************/
pair<TH2F*,TH2F*> Bin2DTree::fillWidths()
/*****************************************************************/
{
        if(!m_gridConstraint)
        {
            cerr<<"ERROR: Bin2DLeaf.fillWidths(): Trying to fill histogram, but the binning is unknown. Define first the gridConstraint.\n";
            exit(1);
        }
        TH2F* hWidthX = (TH2F*)m_gridConstraint->Clone("widthXFromTree");
        TH2F* hWidthY = (TH2F*)m_gridConstraint->Clone("widthYFromTree");
        int nbinsx = hWidthX->GetNbinsX();
        int nbinsy = hWidthX->GetNbinsY();
        map<Bin2DLeaf*, vector< pair<int,int> > > binsInLeaf;
        for(int bx=1;bx<nbinsx+1;bx++) 
        {
            for(int by=1;by<nbinsy+1;by++)
            {
                double x = hWidthX->GetXaxis()->GetBinCenter(bx);
                double y = hWidthX->GetYaxis()->GetBinCenter(by);
                Bin2DLeaf* leaf = getLeaf(x, y);
                vector<Bin2DLeaf*> neighborLeaves = findNeighborLeaves(leaf);
                neighborLeaves.push_back(leaf);
                vector<Bin2DLeaf*>::iterator itLeaf = neighborLeaves.begin();
                vector<Bin2DLeaf*>::iterator itELeaf = neighborLeaves.end();
                double sumw = 0.;
                double sumwx = 0.;
                double sumwy = 0.;
                for(;itLeaf!=itELeaf;++itLeaf)
                {
                    double xi = (*itLeaf)->getXCenter();
                    double yi = (*itLeaf)->getYCenter();
                    double dx = xi-x;
                    double dy = yi-y;
                    double wxi = (*itLeaf)->getXWidth();
                    double wyi = (*itLeaf)->getYWidth();
                    if(dx<0.05*wxi) dx = 0.05*wxi;
                    if(dy<0.05*wyi) dy = 0.05*wyi;
                    double dr2 = dx*dx+dy*dy;
                    sumw += 1./dr2;
                    sumwx += wxi/dr2;
                    sumwy += wyi/dr2;
                }
                double widthx = sumwx/sumw;
                double widthy = sumwy/sumw;
                hWidthX->SetBinContent(bx,by,widthx);
                hWidthY->SetBinContent(bx,by,widthy);
                hWidthX->SetBinError(bx,by,0.);
                hWidthY->SetBinError(bx,by,0.);
            }
        }
        return make_pair(hWidthX, hWidthY);
}





