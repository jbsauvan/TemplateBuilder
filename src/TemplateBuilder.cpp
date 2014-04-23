
#include "TemplateBuilder.h"
#include "BinTree.h"
#include "GaussKernelSmoother.h"
#include "Smoother1D.h"

#include "TH2F.h"
#include "TH3F.h"
#include "TGraph.h"

#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;


/*****************************************************************/
TemplateBuilder::~TemplateBuilder()
/*****************************************************************/
{
    map<string, Template*>::iterator it = m_templates.begin();
    map<string, Template*>::iterator itE = m_templates.end();
    for(;it!=itE;++it)
    {
        delete it->second;
    }
}

/*****************************************************************/
void TemplateBuilder::addTemplate(const std::string& name)
/*****************************************************************/
{
    m_templates[name] = new Template();
    m_templates[name]->setName(name);
}

/*****************************************************************/
void TemplateBuilder::addTemplate(Template* tmp)
/*****************************************************************/
{
    if(m_templates.find(tmp->getName())!=m_templates.end())
    {
        cout<<"[WARN] TemplateBuilder::addTemplate(): template name '"<<tmp->getName()<<"' exists already. The template definition will be overwritten.\n";
    }
    m_templates[tmp->getName()] = tmp;
}


/*****************************************************************/
Template* TemplateBuilder::getTemplate(const std::string& name)
/*****************************************************************/
{
    map<string,Template*>::iterator it = m_templates.find(name);   
    if(it==m_templates.end())
    {
        stringstream error;
        error << "TemplateBuilder::getTemplate(): template '"<<name<<"' doesnt' exist";
        throw runtime_error(error.str());
    }
    return it->second;
}


/*****************************************************************/
void TemplateBuilder::fillTemplates()
/*****************************************************************/
{
    map<string, Template*>::iterator tmpIt = m_templates.begin();
    map<string, Template*>::iterator tmpItE = m_templates.end();
    for(;tmpIt!=tmpItE;++tmpIt)
    {
        Template* tmp = tmpIt->second;
        if(tmp->getOrigin()!=Template::Origin::FILES) continue;

        if(tmp->getBinningType()==Template::BinningType::FIXED)
        {
            cout<< "[INFO] Building "<<tmp->numberOfDimensions()<<"D template '"<<tmp->getName()<<"' with fixed size binning\n";
            int overflows = 0;
            if(tmp->numberOfDimensions()==2)
            {
                TH2F* histo = dynamic_cast<TH2F*>(tmp->getTemplate());
                TH2F* histoRaw = dynamic_cast<TH2F*>(tmp->getRawTemplate());
                for(unsigned int e=0;e<tmp->entries().size();e++)
                {
                    int bin = histo->Fill(tmp->entries()[e][0],tmp->entries()[e][1],tmp->weights()[e]);
                    histoRaw->Fill(tmp->entries()[e][0],tmp->entries()[e][1],tmp->weights()[e]);
                    if(bin!=-1)
                    {
                        tmp->getRaw1DTemplate(0)->Fill(tmp->entries()[e][0], tmp->weights()[e]);
                        tmp->getRaw1DTemplate(1)->Fill(tmp->entries()[e][1], tmp->weights()[e]);
                    }
                    else
                    {
                        overflows++;
                    }
                }
            }
            else if(tmp->numberOfDimensions()==3)
            {
                TH3F* histo = dynamic_cast<TH3F*>(tmp->getTemplate());
                TH3F* histoRaw = dynamic_cast<TH3F*>(tmp->getRawTemplate());
                for(unsigned int e=0;e<tmp->entries().size();e++)
                {
                    int bin = histo->Fill(tmp->entries()[e][0],tmp->entries()[e][1],tmp->entries()[e][2],tmp->weights()[e]);
                    histoRaw->Fill(tmp->entries()[e][0],tmp->entries()[e][1],tmp->entries()[e][2],tmp->weights()[e]);
                    if(bin!=-1)
                    {
                        tmp->getRaw1DTemplate(0)->Fill(tmp->entries()[e][0], tmp->weights()[e]);
                        tmp->getRaw1DTemplate(1)->Fill(tmp->entries()[e][1], tmp->weights()[e]);
                        tmp->getRaw1DTemplate(2)->Fill(tmp->entries()[e][2], tmp->weights()[e]);
                    }
                    else
                    {
                        overflows++;
                    }
                }
            }
            if(overflows>0)
            {
                cout<<"[WARN]   "<<overflows<<" events in under/overflow bins\n";
            }
        }
        else if(tmp->getBinningType()==Template::BinningType::ADAPTIVE)
        {
            cout<< "[INFO] Deriving adaptive binning for "<<tmp->numberOfDimensions()<<"D template '"<<tmp->getName()<<"'\n";
            if(tmp->numberOfDimensions()==2)
            {
                TH2F* histoRaw = dynamic_cast<TH2F*>(tmp->getRawTemplate());
                for(unsigned int e=0;e<tmp->entries().size();e++)
                {
                    histoRaw->Fill(tmp->entries()[e][0],tmp->entries()[e][1],tmp->weights()[e]);
                    tmp->getRaw1DTemplate(0)->Fill(tmp->entries()[e][0], tmp->weights()[e]);
                    tmp->getRaw1DTemplate(1)->Fill(tmp->entries()[e][1], tmp->weights()[e]);
                }
            }
            else if(tmp->numberOfDimensions()==3)
            {
                TH3F* histoRaw = dynamic_cast<TH3F*>(tmp->getRawTemplate());
                for(unsigned int e=0;e<tmp->entries().size();e++)
                {
                    histoRaw->Fill(tmp->entries()[e][0],tmp->entries()[e][1],tmp->entries()[e][2],tmp->weights()[e]);
                    tmp->getRaw1DTemplate(0)->Fill(tmp->entries()[e][0], tmp->weights()[e]);
                    tmp->getRaw1DTemplate(1)->Fill(tmp->entries()[e][1], tmp->weights()[e]);
                    tmp->getRaw1DTemplate(2)->Fill(tmp->entries()[e][2], tmp->weights()[e]);
                }
            }
            BinTree bintree(tmp->getMinMax(), tmp->entries(), tmp->weights());
            bintree.setMinLeafEntries(tmp->getEntriesPerBin());
            TH1* gridConstraint = (TH1*)tmp->getTemplate()->Clone("gridConstraint");
            bintree.setGridConstraint(gridConstraint);
            bintree.build();
            cout<<"[INFO]   Number of bins = "<<bintree.getNLeaves()<<"\n";
            cout<<"[INFO]   Smallest bin widths: wx="<<bintree.getMinBinWidth(0)<<", wy="<<bintree.getMinBinWidth(1);
            if(tmp->numberOfDimensions()==3)
            {
                cout<<", wz="<<bintree.getMinBinWidth(2)<<"\n";
            }
            else
            {
                cout<<"\n";
            }
            TH1* histo = dynamic_cast<TH1*>(bintree.fillHistogram());
            tmp->setTemplate(histo);
            cout<< "[INFO] Computing width maps from adaptive binning for template '"<<tmp->getName()<<"'\n";
            vector<TH1*> widths = bintree.fillWidths();
            tmp->setWidths(widths);
            gridConstraint->Delete();
        }
        // make control plot
        tmp->makeProjectionControlPlot("afterFill");
    }
}


/*****************************************************************/
void TemplateBuilder::postProcessing(Template::Origin origin)
/*****************************************************************/
{
    map<string, Template*>::iterator tmpIt = m_templates.begin();
    map<string, Template*>::iterator tmpItE = m_templates.end();
    for(;tmpIt!=tmpItE;++tmpIt)
    {
        Template* tmp = tmpIt->second;
        if(tmp->getOrigin()!=origin) continue;

        double sumOfweightsBefore = tmp->getTemplate()->GetSumOfWeights();
        double scaleFactor = 1.;

        vector<PostProcessing>::iterator it = tmp->postProcessingBegin();
        vector<PostProcessing>::iterator itE = tmp->postProcessingEnd();
        for(;it!=itE;++it)
        {
            switch(it->type())
            {
                case PostProcessing::Type::SMOOTH:
                    {
                        string kernel = it->getParameter<string>("kernel");
                        if(kernel=="k5b")
                        {
                            if(tmp->numberOfDimensions()!=2)
                            {
                                stringstream error;
                                error << "TemplateBuilder::postProcessing(): ('"<<tmp->getName()<<"') Can only apply k5b smoothing for 2D templates\n";
                                throw runtime_error(error.str());
                            }
                            cout<<"[INFO] Smoothing template '"<<tmp->getName()<<"' with k5b kernel\n";
                            tmp->getTemplate()->Smooth(1, "k5b");
                        }
                        else if(kernel=="adaptive")
                        {
                            cout<<"[INFO] Smoothing template '"<<tmp->getName()<<"' with variable Gaussian kernel\n";
                            // First derive adaptive binning if not already done previously
                            // This is needed to define kernel widths
                            if(tmp->getBinningType()!=Template::BinningType::ADAPTIVE)
                            {
                                vector< pair<double,double> > minmax = tmp->getMinMax();
                                BinTree bintree(minmax, tmp->entries(), tmp->weights());
                                unsigned int entriesPerBin = it->getParameter<unsigned int>("entriesperbin");
                                bintree.setMinLeafEntries(entriesPerBin);
                                cout<< "[INFO]   First deriving "<<tmp->numberOfDimensions()<<"D adaptive binning\n";
                                bintree.build();
                                cout<<"[INFO]   Number of bins = "<<bintree.getNLeaves()<<"\n";
                                cout<<"[INFO]   Smallest bin widths: wx="<<bintree.getMinBinWidth(0)<<", wy="<<bintree.getMinBinWidth(1);
                                if(tmp->numberOfDimensions()==3)
                                {
                                    cout<<", wz="<<bintree.getMinBinWidth(2)<<"\n";
                                }
                                else
                                {
                                    cout<<"\n";
                                }
                                TH1* widthTemplate = (TH1*)tmp->getTemplate()->Clone("widthTemplate");
                                cout<< "[INFO]   Computing width maps from adaptive binning\n";
                                vector<TH1*> widths = bintree.fillWidths(widthTemplate);
                                tmp->setWidths(widths);
                                widthTemplate->Delete();
                                cout<< "[INFO]   Applying smoothing based on the width map\n";
                            }
                            GaussKernelSmoother smoother(tmp->numberOfDimensions());
                            smoother.setWidths(tmp->getWidths());
                            TH1* histoSmooth = smoother.smooth(tmp->getTemplate());
                            tmp->setTemplate(histoSmooth);
                        }
                        else
                        {
                            stringstream error;
                            error << "TemplateBuilder::postProcessing(): ('"<<tmp->getName()<<"') Unknown smoothing kernel '"<<kernel<<"'\n";
                            throw runtime_error(error.str());
                        }
                        tmp->makeProjectionControlPlot("afterSmooth");
                        //tmp->makeResidualsControlPlot("afterSmooth");
                        //tmp->makeResidualsControlPlot("afterSmooth", 2);
                        //tmp->makeResidualsControlPlot("afterSmooth", 5);
                        //tmp->makeResidualsControlPlot("afterSmooth", 10);
                        break;
                    }
                case PostProcessing::Type::MIRROR:
                    {
                        bool antiMirror = it->getParameter<bool>("antisymmetric");
                        unsigned int axis = it->getParameter<unsigned int>("axis");
                        if((unsigned int)axis>=tmp->numberOfDimensions())
                        {
                            stringstream error;
                            error << "TemplateBuilder::postProcessing(): Mirroring "<<tmp->numberOfDimensions()<<"D template '"<<tmp->getName()<<"' along axis "<<axis<<" is not possible (axis numbering starts at 0)\n";
                            throw runtime_error(error.str());
                        }
                        cout<<"[INFO] "<<(antiMirror ? "Anti-mirroring" : "Mirroring")<<" template '"<<tmp->getName()<<"' along axis "<<axis<<"\n";
                        if(tmp->numberOfDimensions()==2)
                        {
                            TH2F* histo = dynamic_cast<TH2F*>(tmp->getTemplate());
                            if(axis==0)
                            {
                                for (int binx=0;binx<histo->GetNbinsX()/2; binx++)
                                {
                                    for (int biny=0;biny<histo->GetNbinsY(); biny++)
                                    {
                                        double avr = (antiMirror ? histo->GetBinContent(binx+1,biny+1) - histo->GetBinContent(histo->GetNbinsX()-binx,biny+1) : histo->GetBinContent(binx+1,biny+1) + histo->GetBinContent(histo->GetNbinsX()-binx,biny+1));
                                        histo->SetBinContent(binx+1, biny+1, avr/2.);
                                        histo->SetBinContent(histo->GetNbinsX()-binx, biny+1, (antiMirror ? -avr/2. : avr/2.));
                                    } 
                                }
                            }
                            else if(axis==1)
                            {
                                for (int binx=0;binx<histo->GetNbinsX(); binx++)
                                {
                                    for (int biny=0;biny<histo->GetNbinsY()/2; biny++)
                                    {
                                        double avr = (antiMirror ? histo->GetBinContent(binx+1,biny+1) - histo->GetBinContent(binx+1,histo->GetNbinsY()-biny) : histo->GetBinContent(binx+1,biny+1) + histo->GetBinContent(binx+1,histo->GetNbinsY()-biny));
                                        histo->SetBinContent(binx+1, biny+1, avr/2.);
                                        histo->SetBinContent(binx+1, histo->GetNbinsY()-biny, (antiMirror ? -avr/2. : avr/2.));
                                    } 
                                }
                            }
                        }
                        else if(tmp->numberOfDimensions()==3)
                        {
                            TH3F* histo = dynamic_cast<TH3F*>(tmp->getTemplate());
                            if(axis==0)
                            {
                                for (int binx=0;binx<histo->GetNbinsX()/2; binx++)
                                {
                                    for (int biny=0;biny<histo->GetNbinsY(); biny++)
                                    {
                                        for (int binz=0;binz<histo->GetNbinsZ(); binz++)
                                        {
                                            double avr = (antiMirror ? histo->GetBinContent(binx+1,biny+1,binz+1) - histo->GetBinContent(histo->GetNbinsX()-binx,biny+1,binz+1) : histo->GetBinContent(binx+1,biny+1,binz+1) + histo->GetBinContent(histo->GetNbinsX()-binx,biny+1,binz+1));
                                            histo->SetBinContent(binx+1, biny+1, binz+1, avr/2.);
                                            histo->SetBinContent(histo->GetNbinsX()-binx, biny+1, biny+1, (antiMirror ? -avr/2. : avr/2.));
                                        }
                                    } 
                                }
                            }
                            else if(axis==1)
                            {
                                for (int binx=0;binx<histo->GetNbinsX(); binx++)
                                {
                                    for (int biny=0;biny<histo->GetNbinsY()/2; biny++)
                                    {
                                        for (int binz=0;binz<histo->GetNbinsZ(); binz++)
                                        {

                                            double avr = (antiMirror ? histo->GetBinContent(binx+1,biny+1,binz+1) - histo->GetBinContent(binx+1,histo->GetNbinsY()-biny,binz+1) : histo->GetBinContent(binx+1,biny+1,binz+1) + histo->GetBinContent(binx+1,histo->GetNbinsY()-biny,binz+1));
                                            histo->SetBinContent(binx+1, biny+1,binz+1, avr/2.);
                                            histo->SetBinContent(binx+1, histo->GetNbinsY()-biny,binz+1, (antiMirror ? -avr/2. : avr/2.));
                                        } 
                                    }
                                }
                            }
                            else if(axis==2)
                            {
                                for (int binx=0;binx<histo->GetNbinsX(); binx++)
                                {
                                    for (int biny=0;biny<histo->GetNbinsY(); biny++)
                                    {
                                        for (int binz=0;binz<histo->GetNbinsZ()/2; binz++)
                                        {
                                            double avr = (antiMirror ? histo->GetBinContent(binx+1,biny+1,binz+1) - histo->GetBinContent(binx+1,biny+1,histo->GetNbinsZ()-binz) : histo->GetBinContent(binx+1,biny+1,binz+1) + histo->GetBinContent(binx+1,biny+1,histo->GetNbinsZ()-binz));
                                            histo->SetBinContent(binx+1, biny+1,binz+1, avr/2.);
                                            histo->SetBinContent(binx+1,biny+1, histo->GetNbinsZ()-binz, (antiMirror ? -avr/2. : avr/2.));
                                        } 
                                    }
                                }
                            }
                        }
                        tmp->makeProjectionControlPlot("afterMirror");
                        //tmp->makeResidualsControlPlot("afterMirror");
                        //tmp->makeResidualsControlPlot("afterMirror", 2);
                        //tmp->makeResidualsControlPlot("afterMirror", 5);
                        //tmp->makeResidualsControlPlot("afterMirror", 10);
                        break;
                    }
                case PostProcessing::Type::FLOOR:
                    {
                        cout<<"[INFO] Flooring template '"<<tmp->getName()<<"'\n";
                        if(tmp->getTemplate()->GetMinimum()>0.)
                        {
                            cout<<"[INFO]   No zero bin. Flooring is not needed.\n";
                            break;
                        }
                        if(tmp->numberOfDimensions()==2)
                        {
                            TH2F* histo = dynamic_cast<TH2F*>(tmp->getTemplate());
                            double floorN = ((histo->Integral())/(histo->GetNbinsX()*histo->GetNbinsY()))*(0.001/100.);
                            for(int binx = 1; binx <= histo->GetNbinsX(); binx++)
                            {
                                for(int biny = 1; biny <= histo->GetNbinsY(); biny++)
                                {
                                    double orig = histo->GetBinContent(binx,biny);
                                    histo->SetBinContent(binx,biny,(orig+floorN));
                                }
                            }
                        }
                        else if(tmp->numberOfDimensions()==3)
                        {
                            TH3F* histo = dynamic_cast<TH3F*>(tmp->getTemplate());
                            double floorN = ((histo->Integral())/(histo->GetNbinsX()*histo->GetNbinsY()*histo->GetNbinsZ()))*(0.001/100.);
                            for(int binx = 1; binx <= histo->GetNbinsX(); binx++)
                            {
                                for(int biny = 1; biny <= histo->GetNbinsY(); biny++)
                                {
                                    for(int binz = 1; binz <= histo->GetNbinsZ(); binz++)
                                    {
                                        double orig = histo->GetBinContent(binx,biny,binz);
                                        histo->SetBinContent(binx,biny,binz,(orig+floorN));
                                    }
                                }
                            }
                        }
                        tmp->makeProjectionControlPlot("afterFloor");
                        //tmp->makeResidualsControlPlot("afterFloor");
                        //tmp->makeResidualsControlPlot("afterFloor", 2);
                        //tmp->makeResidualsControlPlot("afterFloor", 5);
                        //tmp->makeResidualsControlPlot("afterFloor", 10);
                        break;
                    }
                case PostProcessing::Type::RESCALE:
                    {
                        double factor = it->getParameter<double>("factor");
                        scaleFactor *= factor;
                        break;
                    }
                case PostProcessing::Type::REWEIGHT:
                    {
                        cout<<"[INFO] Reweighting template '"<<tmp->getName()<<"'\n";
                        applyReweighting(tmp,*it);
                        tmp->makeProjectionControlPlot("afterReweight");
                        //tmp->makeResidualsControlPlot("afterReweight");
                        //tmp->makeResidualsControlPlot("afterReweight", 2);
                        //tmp->makeResidualsControlPlot("afterReweight", 5);
                        //tmp->makeResidualsControlPlot("afterReweight", 10);
                        break;
                    }
                default:
                    break;
            }
            if(it->type()!=PostProcessing::Type::RESCALE)
            {
                double sumOfWeightsAfter = tmp->getTemplate()->GetSumOfWeights();
                cout<<"[INFO]   Sum of weights after/before = "<<sumOfWeightsAfter<<" / "<<sumOfweightsBefore<<" = "<<sumOfWeightsAfter/sumOfweightsBefore<<"\n";
                sumOfweightsBefore = sumOfWeightsAfter;
            }
        }
        // normalize
        if(origin==Template::Origin::FILES)
        {
            double targetSumOfWeights = 1.;
            if(tmp->conserveSumOfWeights())
            {
                targetSumOfWeights = tmp->originalSumOfWeights();
                cout<<"[INFO] Normalizing template '"<<tmp->getName()<<"' to the original sum of weights = "<<targetSumOfWeights<<"\n";
            }
            else
            {
                cout<<"[INFO] Normalizing template '"<<tmp->getName()<<"' to 1\n";
            }
            double sumOfWeights = tmp->getTemplate()->GetSumOfWeights();
            tmp->getTemplate()->Scale(targetSumOfWeights/sumOfWeights);
            sumOfWeights = tmp->getRawTemplate()->GetSumOfWeights();
            tmp->getRawTemplate()->Scale(targetSumOfWeights/sumOfWeights);
            for(unsigned int axis=0;axis<tmp->numberOfDimensions();axis++)
            {
                sumOfWeights = tmp->getRaw1DTemplate(axis)->GetSumOfWeights();
                tmp->getRaw1DTemplate(axis)->Scale(targetSumOfWeights/sumOfWeights);
            }
        }
        //double scaleFactor = tmp->getRescaling();
        if(scaleFactor!=1.)
        {
            cout<<"[INFO] Rescaling template '"<<tmp->getName()<<"' with factor "<<scaleFactor<<"\n";
            tmp->getTemplate()->Scale(scaleFactor);
            tmp->getRawTemplate()->Scale(scaleFactor);
            for(unsigned int axis=0;axis<tmp->numberOfDimensions();axis++)
            {
                tmp->getRaw1DTemplate(axis)->Scale(scaleFactor);
            }
        }
        tmp->makeProjectionControlPlot("afterNormalization");
        //tmp->makeResidualsControlPlot("afterNormalization");
        //tmp->makeResidualsControlPlot("afterNormalization", 2);
        //tmp->makeResidualsControlPlot("afterNormalization", 5);
        //tmp->makeResidualsControlPlot("afterNormalization", 10);
    }

}

/*****************************************************************/
void TemplateBuilder::buildTemplatesFromTemplates()
/*****************************************************************/
{
    map<string, Template*>::iterator tmpIt = m_templates.begin();
    map<string, Template*>::iterator tmpItE = m_templates.end();
    for(;tmpIt!=tmpItE;++tmpIt)
    {
        Template* tmp = tmpIt->second;
        if(tmp->getOrigin()!=Template::Origin::TEMPLATES) continue;

        cout<<"[INFO] Building template '"<<tmp->getName()<<"' from previously filled templates\n";

        vector<pair<string, double> >::const_iterator it = tmp->inputTemplatesBegin();
        vector<pair<string, double> >::const_iterator itE = tmp->inputTemplatesEnd();
        for(;it!=itE;++it)
        {
            string inTmpName = it->first;
            double factor = it->second;
            map<string, Template*>::const_iterator inTmpIt = m_templates.find(inTmpName);
            if(inTmpIt==m_templates.end())
            {
                stringstream error;
                error << "TemplateBuilder::buildTemplatesFromTemplates(): Cannot find template '"<<inTmpName<<"'\n";
                throw runtime_error(error.str());
            }
            Template* inTmp = inTmpIt->second;
            if(!tmp->getTemplate())
            {
                cout<<"[INFO]   + ("<<factor<<") x "<<inTmp->getName()<<"\n";
                tmp->setTemplate(inTmp->getTemplate());
                tmp->setRawTemplate(inTmp->getRawTemplate());
                tmp->getTemplate()->Scale(factor);
                tmp->getRawTemplate()->Scale(factor);
                tmp->setWidths(inTmp->getWidths());
                tmp->setRaw1DTemplates(inTmp->getRaw1DTemplates());
                tmp->setRaw2DTemplates(inTmp->getRaw2DTemplates());
                for(unsigned int axis=0;axis<inTmp->numberOfDimensions();axis++)
                {
                    tmp->getRaw1DTemplate(axis)->Scale(factor);
                }
                for(unsigned int i=0;i<tmp->getRaw2DTemplates().size();i++)
                {
                    tmp->getRaw2DTemplate(i)->Scale(factor);
                }
                vector<string> vars;
                for(unsigned int v=0;v<inTmp->numberOfDimensions();v++)
                {
                    vars.push_back(inTmp->getVariable(v));
                }
                tmp->setVariables(vars);
            }
            else
            {
                // TODO: take the averaged width for each bin instead of the width of the first template
                cout<<"[INFO]   + ("<<factor<<") x "<<inTmp->getName()<<"\n";
                // FIXME: check this for N dimensions
                //int nbinsx1 = tmp->getTemplate()->GetNbinsX();
                //int nbinsy1 = tmp->getTemplate()->GetNbinsY();
                //double minx1 = tmp->getTemplate()->GetXaxis()->GetBinLowEdge(0);
                //double maxx1 = tmp->getTemplate()->GetXaxis()->GetBinUpEdge(nbinsx1);
                //double miny1 = tmp->getTemplate()->GetYaxis()->GetBinLowEdge(0);
                //double maxy1 = tmp->getTemplate()->GetYaxis()->GetBinUpEdge(nbinsy1);
                //int nbinsx2 = inTmp->getTemplate()->GetNbinsX();
                //int nbinsy2 = inTmp->getTemplate()->GetNbinsY();
                //double minx2 = inTmp->getTemplate()->GetXaxis()->GetBinLowEdge(0);
                //double maxx2 = inTmp->getTemplate()->GetXaxis()->GetBinUpEdge(nbinsx2);
                //double miny2 = inTmp->getTemplate()->GetYaxis()->GetBinLowEdge(0);
                //double maxy2 = inTmp->getTemplate()->GetYaxis()->GetBinUpEdge(nbinsy2);
                //if(nbinsx1!=nbinsx2 || nbinsy1!=nbinsy2 ||
                //        minx1!=minx2 || maxx1!=maxx2 ||
                //        miny1!=miny2 || maxy1!=maxy2
                //        )
                //{
                //    stringstream error;
                //    error << "TemplateBuilder::buildTemplatesFromTemplates(): Trying to add templates with different binning ("<<tmp->getName()<<" & "<<inTmp->getName()<<")\n";
                //    throw runtime_error(error.str());
                //}
                tmp->getTemplate()->Add(inTmp->getTemplate(), factor);
                tmp->getRawTemplate()->Add(inTmp->getRawTemplate(), factor);
                for(unsigned int axis=0;axis<inTmp->numberOfDimensions();axis++)
                {
                    tmp->getRaw1DTemplate(axis)->Add(inTmp->getRaw1DTemplate(axis), factor);
                }
                for(unsigned int i=0;i<tmp->getRaw2DTemplates().size();i++)
                {
                    tmp->getRaw2DTemplate(i)->Scale(factor);
                }
            }
        }
    }


    /// postprocessing
    //tmpIt = m_templates.begin();
    //tmpItE = m_templates.end();
    //for(;tmpIt!=tmpItE;++tmpIt)
    //{
    //    Template* tmp = tmpIt->second;
    //    if(tmp->getOrigin()!=Template::Origin::TEMPLATES) continue;

    //    double sumOfweightsBefore = tmp->getTemplate()->GetSumOfWeights();
    //    double scaleFactor = 1.;

    //    vector<PostProcessing>::iterator it = tmp->postProcessingBegin();
    //    vector<PostProcessing>::iterator itE = tmp->postProcessingEnd();
    //    for(;it!=itE;++it)
    //    {
    //        switch(it->type())
    //        {
    //            case PostProcessing::Type::SMOOTH:
    //                {
    //                    string kernel = "";
    //                    it->getParameter("kernel", kernel);
    //                    if(kernel=="k5b")
    //                    {
    //                        if(tmp->numberOfDimensions()!=2)
    //                        {
    //                            stringstream error;
    //                            error << "TemplateBuilder::postProcessing(): ('"<<tmp->getName()<<"') Can only apply k5b smoothing for 2D templates\n";
    //                            throw runtime_error(error.str());
    //                        }
    //                        cout<<"[INFO] Smoothing template '"<<tmp->getName()<<"' with k5b kernel\n";
    //                        tmp->getTemplate()->Smooth(1, "k5b");
    //                    }
    //                    else if(kernel=="adaptive")
    //                    {
    //                        cout<<"[INFO] Smoothing template '"<<tmp->getName()<<"' with variable Gaussian kernel\n";
    //                        // First derive adaptive binning if not already done previously
    //                        // This is needed to define kernel widths
    //                        if(tmp->getBinningType()!=Template::BinningType::ADAPTIVE)
    //                        {
    //                            vector< pair<double,double> > minmax = tmp->getMinMax();
    //                            BinTree bintree(minmax, tmp->entries(), tmp->weights());
    //                            int entriesPerBin = 200;
    //                            it->getParameter("entriesperbin", entriesPerBin);
    //                            bintree.setMinLeafEntries(entriesPerBin);
    //                            cout<< "[INFO]   First deriving "<<tmp->numberOfDimensions()<<"D adaptive binning\n";
    //                            bintree.build();
    //                            cout<<"[INFO]   Number of bins = "<<bintree.getNLeaves()<<"\n";
    //                            cout<<"[INFO]   Smallest bin widths: wx="<<bintree.getMinBinWidth(0)<<", wy="<<bintree.getMinBinWidth(1);
    //                            if(tmp->numberOfDimensions()==3)
    //                            {
    //                                cout<<", wz="<<bintree.getMinBinWidth(2)<<"\n";
    //                            }
    //                            else
    //                            {
    //                                cout<<"\n";
    //                            }
    //                            TH1* widthTemplate = (TH1*)tmp->getTemplate()->Clone("widthTemplate");
    //                            cout<< "[INFO]   Computing width maps from adaptive binning\n";
    //                            vector<TH1*> widths = bintree.fillWidths(widthTemplate);
    //                            tmp->setWidths(widths);
    //                            widthTemplate->Delete();
    //                            cout<< "[INFO]   Applying smoothing based on the width map\n";
    //                        }
    //                        GaussKernelSmoother smoother(tmp->numberOfDimensions());
    //                        smoother.setWidths(tmp->getWidths());
    //                        TH1* histoSmooth = smoother.smooth(tmp->getTemplate());
    //                        tmp->setTemplate(histoSmooth);
    //                    }
    //                    else
    //                    {
    //                        stringstream error;
    //                        error << "TemplateBuilder::postProcessing(): ('"<<tmp->getName()<<"') Unknown smoothing kernel '"<<kernel<<"'\n";
    //                        throw runtime_error(error.str());
    //                    }
    //                    break;
    //                }
    //            case PostProcessing::Type::MIRROR:
    //                {
    //                    bool antiMirror = false;
    //                    it->getParameter("antisymmetric", antiMirror);
    //                    int axis = 1;
    //                    it->getParameter("axis", axis);
    //                    if((unsigned int)axis>=tmp->numberOfDimensions())
    //                    {
    //                        stringstream error;
    //                        error << "TemplateBuilder::postProcessing(): Mirroring "<<tmp->numberOfDimensions()<<"D template '"<<tmp->getName()<<"' along axis "<<axis<<" is not possible (axis numbering starts at 0)\n";
    //                        throw runtime_error(error.str());
    //                    }
    //                    cout<<"[INFO] "<<(antiMirror ? "Anti-mirroring" : "Mirroring")<<" template '"<<tmp->getName()<<"' along axis "<<axis<<"\n";
    //                    if(tmp->numberOfDimensions()==2)
    //                    {
    //                        TH2F* histo = dynamic_cast<TH2F*>(tmp->getTemplate());
    //                        if(axis==0)
    //                        {
    //                            for (int binx=0;binx<histo->GetNbinsX()/2; binx++)
    //                            {
    //                                for (int biny=0;biny<histo->GetNbinsY(); biny++)
    //                                {
    //                                    double avr = (antiMirror ? histo->GetBinContent(binx+1,biny+1) - histo->GetBinContent(histo->GetNbinsX()-binx,biny+1) : histo->GetBinContent(binx+1,biny+1) + histo->GetBinContent(histo->GetNbinsX()-binx,biny+1));
    //                                    histo->SetBinContent(binx+1, biny+1, avr/2.);
    //                                    histo->SetBinContent(histo->GetNbinsX()-binx, biny+1, (antiMirror ? -avr/2. : avr/2.));
    //                                } 
    //                            }
    //                        }
    //                        else if(axis==1)
    //                        {
    //                            for (int binx=0;binx<histo->GetNbinsX(); binx++)
    //                            {
    //                                for (int biny=0;biny<histo->GetNbinsY()/2; biny++)
    //                                {
    //                                    double avr = (antiMirror ? histo->GetBinContent(binx+1,biny+1) - histo->GetBinContent(binx+1,histo->GetNbinsY()-biny) : histo->GetBinContent(binx+1,biny+1) + histo->GetBinContent(binx+1,histo->GetNbinsY()-biny));
    //                                    histo->SetBinContent(binx+1, biny+1, avr/2.);
    //                                    histo->SetBinContent(binx+1, histo->GetNbinsY()-biny, (antiMirror ? -avr/2. : avr/2.));
    //                                } 
    //                            }
    //                        }
    //                    }
    //                    else if(tmp->numberOfDimensions()==3)
    //                    {
    //                        TH3F* histo = dynamic_cast<TH3F*>(tmp->getTemplate());
    //                        if(axis==0)
    //                        {
    //                            for (int binx=0;binx<histo->GetNbinsX()/2; binx++)
    //                            {
    //                                for (int biny=0;biny<histo->GetNbinsY(); biny++)
    //                                {
    //                                    for (int binz=0;binz<histo->GetNbinsZ(); binz++)
    //                                    {
    //                                        double avr = (antiMirror ? histo->GetBinContent(binx+1,biny+1,binz+1) - histo->GetBinContent(histo->GetNbinsX()-binx,biny+1,binz+1) : histo->GetBinContent(binx+1,biny+1,binz+1) + histo->GetBinContent(histo->GetNbinsX()-binx,biny+1,binz+1));
    //                                        histo->SetBinContent(binx+1, biny+1, binz+1, avr/2.);
    //                                        histo->SetBinContent(histo->GetNbinsX()-binx, biny+1, biny+1, (antiMirror ? -avr/2. : avr/2.));
    //                                    }
    //                                } 
    //                            }
    //                        }
    //                        else if(axis==1)
    //                        {
    //                            for (int binx=0;binx<histo->GetNbinsX(); binx++)
    //                            {
    //                                for (int biny=0;biny<histo->GetNbinsY()/2; biny++)
    //                                {
    //                                    for (int binz=0;binz<histo->GetNbinsZ(); binz++)
    //                                    {

    //                                        double avr = (antiMirror ? histo->GetBinContent(binx+1,biny+1,binz+1) - histo->GetBinContent(binx+1,histo->GetNbinsY()-biny,binz+1) : histo->GetBinContent(binx+1,biny+1,binz+1) + histo->GetBinContent(binx+1,histo->GetNbinsY()-biny,binz+1));
    //                                        histo->SetBinContent(binx+1, biny+1,binz+1, avr/2.);
    //                                        histo->SetBinContent(binx+1, histo->GetNbinsY()-biny,binz+1, (antiMirror ? -avr/2. : avr/2.));
    //                                    } 
    //                                }
    //                            }
    //                        }
    //                        else if(axis==2)
    //                        {
    //                            for (int binx=0;binx<histo->GetNbinsX(); binx++)
    //                            {
    //                                for (int biny=0;biny<histo->GetNbinsY(); biny++)
    //                                {
    //                                    for (int binz=0;binz<histo->GetNbinsZ()/2; binz++)
    //                                    {
    //                                        double avr = (antiMirror ? histo->GetBinContent(binx+1,biny+1,binz+1) - histo->GetBinContent(binx+1,biny+1,histo->GetNbinsZ()-binz) : histo->GetBinContent(binx+1,biny+1,binz+1) + histo->GetBinContent(binx+1,biny+1,histo->GetNbinsZ()-binz));
    //                                        histo->SetBinContent(binx+1, biny+1,binz+1, avr/2.);
    //                                        histo->SetBinContent(binx+1,biny+1, histo->GetNbinsZ()-binz, (antiMirror ? -avr/2. : avr/2.));
    //                                    } 
    //                                }
    //                            }
    //                        }
    //                    }
    //                    break;
    //                }
    //            case PostProcessing::Type::FLOOR:
    //                {
    //                    cout<<"[INFO] Flooring template '"<<tmp->getName()<<"'\n";
    //                    if(tmp->getTemplate()->GetMinimum()>0.)
    //                    {
    //                        cout<<"[INFO]   No zero bin. Flooring is not needed.\n";
    //                        break;
    //                    }
    //                    if(tmp->numberOfDimensions()==2)
    //                    {
    //                        TH2F* histo = dynamic_cast<TH2F*>(tmp->getTemplate());
    //                        double floorN = ((histo->Integral())/(histo->GetNbinsX()*histo->GetNbinsY()))*(0.001/100.);
    //                        for(int binx = 1; binx <= histo->GetNbinsX(); binx++)
    //                        {
    //                            for(int biny = 1; biny <= histo->GetNbinsY(); biny++)
    //                            {
    //                                double orig = histo->GetBinContent(binx,biny);
    //                                histo->SetBinContent(binx,biny,(orig+floorN));
    //                            }
    //                        }
    //                    }
    //                    else if(tmp->numberOfDimensions()==3)
    //                    {
    //                        TH3F* histo = dynamic_cast<TH3F*>(tmp->getTemplate());
    //                        double floorN = ((histo->Integral())/(histo->GetNbinsX()*histo->GetNbinsY()*histo->GetNbinsZ()))*(0.001/100.);
    //                        for(int binx = 1; binx <= histo->GetNbinsX(); binx++)
    //                        {
    //                            for(int biny = 1; biny <= histo->GetNbinsY(); biny++)
    //                            {
    //                                for(int binz = 1; binz <= histo->GetNbinsZ(); binz++)
    //                                {
    //                                    double orig = histo->GetBinContent(binx,biny,binz);
    //                                    histo->SetBinContent(binx,biny,binz,(orig+floorN));
    //                                }
    //                            }
    //                        }
    //                    }
    //                    break;
    //                }
    //            case PostProcessing::Type::RESCALE:
    //                {
    //                    double factor = 1.;
    //                    it->getParameter("factor", factor);
    //                    scaleFactor *= factor;
    //                    break;
    //                }
    //            default:
    //                break;
    //        }
    //        double sumOfWeightsAfter = tmp->getTemplate()->GetSumOfWeights();
    //        cout<<"[INFO]   Sum of weights after/before = "<<sumOfWeightsAfter<<" / "<<sumOfweightsBefore<<" = "<<sumOfWeightsAfter/sumOfweightsBefore<<"\n";
    //        sumOfweightsBefore = sumOfWeightsAfter;
    //    }
    //    // normalize
    //    //cout<<"[INFO] Normalizing template '"<<tmp->getName()<<"' to 1\n";
    //    //double sumOfWeights = tmp->getTemplate()->GetSumOfWeights();
    //    //tmp->getTemplate()->Scale(1./sumOfWeights);
    //    //double scaleFactor = tmp->getRescaling();
    //    if(scaleFactor!=1.)
    //    {
    //        tmp->getTemplate()->Scale(scaleFactor);
    //        cout<<"[INFO] Rescaling template '"<<tmp->getName()<<"' with factor "<<scaleFactor<<"\n";
    //    }
    //}
}


/*****************************************************************/
void TemplateBuilder::applyReweighting(Template* tmp, const PostProcessing& pp)
/*****************************************************************/
{
    vector<unsigned int> axes = pp.getParameter< vector<unsigned int> >("axes");
    vector< vector<double> > rebinning = pp.getParameter< vector< vector<double> > >("rebinning");
    vector<unsigned int>::const_iterator it = axes.begin();
    vector<unsigned int>::const_iterator itE = axes.end();
    for(;it!=itE;++it)
    {
        unsigned int axis = *it;
        cout<<"[INFO]   Reweighting along axis "<<axis<<"\n";
        TH1D* refHisto = tmp->getRaw1DTemplate(axis);
        TH1D* projTmp = tmp->getProjected1DTemplate(axis);
        // rebin
        bool rebin = false;
        if(rebinning.size()>axis && rebinning[axis].size()>2)
        {
            rebin = true;
            unsigned int nbins = rebinning[axis].size()-1;
            double* bins = new double[nbins+1];
            for(unsigned int i=0;i<nbins+1;i++)
            {
                bins[i] = rebinning[axis][i];
            }
            TH1D* refHistoRebin = dynamic_cast<TH1D*>(refHisto->Rebin(nbins, "refRebinHisto", bins));
            // rescale according to bin widths
            double oldBinWidth = refHisto->GetXaxis()->GetBinWidth(1);
            for(int b=1;b<refHistoRebin->GetNbinsX()+1;b++)
            {
                double content     = refHistoRebin->GetBinContent(b);
                double error       = refHistoRebin->GetBinError(b);
                double newBinWidth = refHistoRebin->GetXaxis()->GetBinWidth(b);
                refHistoRebin->SetBinContent(b,content/newBinWidth*oldBinWidth);
                refHistoRebin->SetBinError(b,error/newBinWidth*oldBinWidth);
            }
            TGraph* refGraph = new TGraph(refHistoRebin);
            TH1D* refHistoSmooth = dynamic_cast<TH1D*>(refHisto->Clone("refHisto"));
            for(int b=1;b<refHisto->GetNbinsX()+1;b++)
            {
                double x = refHisto->GetXaxis()->GetBinCenter(b);
                double content = max(0.,refGraph->Eval(x, 0, "S"));
                refHistoSmooth->SetBinContent(b, content);
            }
            refHisto = refHistoSmooth;
            delete refHistoRebin;
            delete refGraph;
        }
        else // No rebinning. Automatic procedure
        {
            Smoother1D smoother;
            refHisto = smoother.smooth(refHisto);
        }
        unsigned int nbins = refHisto->GetNbinsX();
        for(unsigned int b=1;b<=nbins;b++)
        {
            double ref = refHisto->GetBinContent(b);
            double old = projTmp->GetBinContent(b);
            double weight = (old!=0. ? ref/old : 1.);
            tmp->reweight1D(axis, b, weight);
        }
        if(rebin)
        {
            delete refHisto;
        }
        delete projTmp;
    }
}
