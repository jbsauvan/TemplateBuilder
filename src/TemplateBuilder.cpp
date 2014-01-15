
#include "TemplateBuilder.h"
#include "BinTree.h"
#include "GaussKernelSmoother.h"

#include "TH2F.h"
#include "TH3F.h"

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
        cout<<"[WARNING] TemplateBuilder::addTemplate(): template name '"<<tmp->getName()<<"' exists already. The template definition will be overwritten.\n";
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
            cout<< "[INFO] Building "<<tmp->numberOfDimensions()<<"D template '"<<tmp->getName()<<"' with standard binning\n";
            if(tmp->numberOfDimensions()==2)
            {
                TH2F* histo = dynamic_cast<TH2F*>(tmpIt->second->getTemplate());
                for(unsigned int e=0;e<tmp->entries().size();e++)
                {
                    histo->Fill(tmp->entries()[e][0],tmp->entries()[e][1],tmp->weights()[e]);
                }
            }
            else if(tmp->numberOfDimensions()==3)
            {
                TH3F* histo = dynamic_cast<TH3F*>(tmpIt->second->getTemplate());
                for(unsigned int e=0;e<tmp->entries().size();e++)
                {
                    histo->Fill(tmp->entries()[e][0],tmp->entries()[e][1],tmp->entries()[e][2],tmp->weights()[e]);
                }
            }
        }
        else if(tmp->getBinningType()==Template::BinningType::ADAPTIVE)
        {
            cout<< "[INFO] Deriving adaptive binning for "<<tmp->numberOfDimensions()<<"D template '"<<tmp->getName()<<"'\n";
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
    }
}


/*****************************************************************/
void TemplateBuilder::postProcessing()
/*****************************************************************/
{
    map<string, Template*>::iterator tmpIt = m_templates.begin();
    map<string, Template*>::iterator tmpItE = m_templates.end();
    for(;tmpIt!=tmpItE;++tmpIt)
    {
        Template* tmp = tmpIt->second;
        if(tmp->getOrigin()!=Template::Origin::FILES) continue;

        double sumOfweightsBefore = tmp->getTemplate()->GetSumOfWeights();

        vector<Template::PostProcessing>::iterator it = tmp->postProcessingBegin();
        vector<Template::PostProcessing>::iterator itE = tmp->postProcessingEnd();
        for(;it!=itE;++it)
        {
            switch(*it)
            {
                case Template::PostProcessing::SMOOTH_K5B:
                    {
                        if(tmp->numberOfDimensions()!=2)
                        {
                            stringstream error;
                            error << "TemplateBuilder::postProcessing(): ('"<<tmp->getName()<<"') Can only apply k5b smoothing for histo with 2D\n";
                            throw runtime_error(error.str());
                        }
                        cout<<"[INFO] Smoothing template '"<<tmp->getName()<<"' with k5b kernel\n";
                        tmp->getTemplate()->Smooth(1, "k5b");
                        break;
                    }
                case Template::PostProcessing::SMOOTH_ADAPTIVE:
                    {
                        cout<<"[INFO] Smoothing template '"<<tmp->getName()<<"' with variable Gaussian kernel\n";
                        // First derive adaptive binning if not already done previously
                        // This is needed to define kernel widths
                        if(tmp->getBinningType()!=Template::BinningType::ADAPTIVE)
                        {
                            vector< pair<double,double> > minmax = tmp->getMinMax();
                            BinTree bintree(minmax, tmp->entries(), tmp->weights());
                            bintree.setMinLeafEntries(tmp->getEntriesPerBin());
                            bintree.build();
                            TH1* widthTemplate = (TH1*)tmp->getTemplate()->Clone("widthTemplate");
                            vector<TH1*> widths = bintree.fillWidths(widthTemplate);

                            tmp->setWidths(widths);
                            widthTemplate->Delete();
                        }
                        GaussKernelSmoother smoother(tmp->numberOfDimensions());
                        smoother.setWidths(tmp->getWidths());
                        TH1* histoSmooth = smoother.smooth(tmp->getTemplate());
                        tmp->setTemplate(histoSmooth);
                        break;
                    }
                case Template::PostProcessing::MIRROR:
                    {
                        cout<<"[INFO] Mirroring template '"<<tmp->getName()<<"'\n";
                        cout<<"[INFO]   !! Mirror is automatically applied on the 2nd dimension for the moment !!\n";
                        // FIXME: mirror can only be done on the second dimension
                        if(tmp->numberOfDimensions()==2)
                        {
                            TH2F* histo = dynamic_cast<TH2F*>(tmp->getTemplate());
                            for (int binx=0;binx<histo->GetNbinsX(); binx++)
                            {
                                for (int biny=0;biny<histo->GetNbinsY()/2; biny++)
                                {
                                    double avr = histo->GetBinContent(binx+1,biny+1) + histo->GetBinContent(binx+1,histo->GetNbinsY()-biny);
                                    histo->SetBinContent(binx+1, biny+1, avr/2.);
                                    histo->SetBinContent(binx+1, histo->GetNbinsY()-biny, avr/2.);
                                } 
                            }
                        }
                        else if(tmp->numberOfDimensions()==3)
                        {
                            TH3F* histo = dynamic_cast<TH3F*>(tmp->getTemplate());
                            for (int binx=0;binx<histo->GetNbinsX(); binx++)
                            {
                                for (int binz=0;binz<histo->GetNbinsZ(); binz++)
                                {
                                    for (int biny=0;biny<histo->GetNbinsY()/2; biny++)
                                    {
                                        double avr = histo->GetBinContent(binx+1,biny+1,binz+1) + histo->GetBinContent(binx+1,histo->GetNbinsY()-biny,binz+1);
                                        histo->SetBinContent(binx+1, biny+1, binz+1, avr/2.);
                                        histo->SetBinContent(binx+1, histo->GetNbinsY()-biny, binz+1, avr/2.);
                                    } 
                                }
                            }
                        }
                        break;
                    }
                case Template::PostProcessing::MIRROR_INV:
                    {
                        cout<<"[INFO] Anti-mirroring template '"<<tmp->getName()<<"'\n";
                        cout<<"[INFO]   !! Anti-mirror is automatically be applied on the 2nd dimension for the moment !!\n";
                        // FIXME: mirror can only be done on the second dimension
                        if(tmp->numberOfDimensions()==2)
                        {
                            TH2F* histo = dynamic_cast<TH2F*>(tmp->getTemplate());
                            for (int binx=0;binx<histo->GetNbinsX(); binx++)
                            {
                                for (int biny=0;biny<histo->GetNbinsY()/2; biny++)
                                {
                                    double avr = histo->GetBinContent(binx+1,biny+1) - histo->GetBinContent(binx+1,histo->GetNbinsY()-biny);
                                    histo->SetBinContent(binx+1, biny+1, avr/2.);
                                    histo->SetBinContent(binx+1, histo->GetNbinsY()-biny, -avr/2.);
                                } 
                            }
                        }
                        else if(tmp->numberOfDimensions()==3)
                        {
                            TH3F* histo = dynamic_cast<TH3F*>(tmp->getTemplate());
                            for (int binx=0;binx<histo->GetNbinsX(); binx++)
                            {
                                for (int binz=0;binz<histo->GetNbinsZ(); binz++)
                                {
                                    for (int biny=0;biny<histo->GetNbinsY()/2; biny++)
                                    {
                                        double avr = histo->GetBinContent(binx+1,biny+1,binz+1) - histo->GetBinContent(binx+1,histo->GetNbinsY()-biny,binz+1);
                                        histo->SetBinContent(binx+1, biny+1, binz+1, avr/2.);
                                        histo->SetBinContent(binx+1, histo->GetNbinsY()-biny, binz+1, -avr/2.);
                                    } 
                                }
                            }
                        }
                        break;
                    }
                case Template::PostProcessing::FLOOR:
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
                        break;
                    }
                default:
                    break;
            }
            double sumOfWeightsAfter = tmp->getTemplate()->GetSumOfWeights();
            cout<<"[INFO]   Sum of weights after/before = "<<sumOfWeightsAfter<<" / "<<sumOfweightsBefore<<" = "<<sumOfWeightsAfter/sumOfweightsBefore<<"\n";
            sumOfweightsBefore = sumOfWeightsAfter;
        }
        // normalize
        cout<<"[INFO] Normalizing template '"<<tmp->getName()<<"' to 1\n";
        double sumOfWeights = tmp->getTemplate()->GetSumOfWeights();
        tmp->getTemplate()->Scale(1./sumOfWeights);
        double scaleFactor = tmp->getRescaling();
        if(scaleFactor!=1.)
        {
            tmp->getTemplate()->Scale(scaleFactor);
            cout<<"[INFO] Rescaling template '"<<tmp->getName()<<"' with factor "<<scaleFactor<<"\n";
        }
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
                cout<<"[INFO] + ("<<factor<<") x "<<inTmp->getName()<<"\n";
                tmp->setTemplate(inTmp->getTemplate());
                tmp->getTemplate()->Scale(factor);
                tmp->setWidths(inTmp->getWidths());
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
                cout<<"[INFO] + ("<<factor<<") x "<<inTmp->getName()<<"\n";
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
            }
        }
    }


    /// postprocessing
    tmpIt = m_templates.begin();
    tmpItE = m_templates.end();
    for(;tmpIt!=tmpItE;++tmpIt)
    {
        Template* tmp = tmpIt->second;
        if(tmp->getOrigin()!=Template::Origin::TEMPLATES) continue;

        double sumOfweightsBefore = tmp->getTemplate()->GetSumOfWeights();

        vector<Template::PostProcessing>::iterator it = tmp->postProcessingBegin();
        vector<Template::PostProcessing>::iterator itE = tmp->postProcessingEnd();
        for(;it!=itE;++it)
        {
            switch(*it)
            {
                case Template::PostProcessing::SMOOTH_K5B:
                    {
                        if(tmp->numberOfDimensions()!=2)
                        {
                            stringstream error;
                            error << "TemplateBuilder::postProcessing(): ('"<<tmp->getName()<<"') Can only apply k5b smoothing for histo with 2D\n";
                            throw runtime_error(error.str());
                        }
                        cout<<"[INFO] Smoothing template '"<<tmp->getName()<<"' with k5b kernel\n";
                        tmp->getTemplate()->Smooth(1, "k5b");
                        break;
                    }
                case Template::PostProcessing::SMOOTH_ADAPTIVE:
                    {
                        cout<<"[INFO] Smoothing template '"<<tmp->getName()<<"' with variable Gaussian kernel\n";
                        // First derive adaptive binning if not already done previously
                        // This is needed to define kernel widths
                        if(tmp->getBinningType()!=Template::BinningType::ADAPTIVE)
                        {
                            vector< pair<double,double> > minmax = tmp->getMinMax();
                            BinTree bintree(minmax, tmp->entries(), tmp->weights());
                            bintree.setMinLeafEntries(tmp->getEntriesPerBin());
                            bintree.build();
                            TH1* widthTemplate = (TH1*)tmp->getTemplate()->Clone("widthTemplate");
                            vector<TH1*> widths = bintree.fillWidths(widthTemplate);
                            tmp->setWidths(widths);
                            widthTemplate->Delete();
                        }
                        GaussKernelSmoother smoother(tmp->numberOfDimensions());
                        smoother.setWidths(tmp->getWidths());
                        TH1* histoSmooth = smoother.smooth(tmp->getTemplate());
                        tmp->setTemplate(histoSmooth);
                        break;
                    }
                case Template::PostProcessing::MIRROR:
                    {
                        cout<<"[INFO] Mirroring template '"<<tmp->getName()<<"'\n";
                        cout<<"[INFO]   !! Mirror is automatically applied on the 2nd dimension for the moment !!\n";
                        // FIXME: mirror can only be done on the second dimension
                        if(tmp->numberOfDimensions()==2)
                        {
                            TH2F* histo = dynamic_cast<TH2F*>(tmp->getTemplate());
                            for (int binx=0;binx<histo->GetNbinsX(); binx++)
                            {
                                for (int biny=0;biny<histo->GetNbinsY()/2; biny++)
                                {
                                    double avr = histo->GetBinContent(binx+1,biny+1) + histo->GetBinContent(binx+1,histo->GetNbinsY()-biny);
                                    histo->SetBinContent(binx+1, biny+1, avr/2.);
                                    histo->SetBinContent(binx+1, histo->GetNbinsY()-biny, avr/2.);
                                } 
                            }
                        }
                        else if(tmp->numberOfDimensions()==3)
                        {
                            TH3F* histo = dynamic_cast<TH3F*>(tmp->getTemplate());
                            for (int binx=0;binx<histo->GetNbinsX(); binx++)
                            {
                                for (int binz=0;binz<histo->GetNbinsZ(); binz++)
                                {
                                    for (int biny=0;biny<histo->GetNbinsY()/2; biny++)
                                    {
                                        double avr = histo->GetBinContent(binx+1,biny+1,binz+1) + histo->GetBinContent(binx+1,histo->GetNbinsY()-biny,binz+1);
                                        histo->SetBinContent(binx+1, biny+1, binz+1, avr/2.);
                                        histo->SetBinContent(binx+1, histo->GetNbinsY()-biny, binz+1, avr/2.);
                                    } 
                                }
                            }
                        }
                        break;
                    }
                case Template::PostProcessing::MIRROR_INV:
                    {
                        cout<<"[INFO] Anti-mirroring template '"<<tmp->getName()<<"'\n";
                        cout<<"[INFO]   !! Anti-mirror is automatically applied on the 2nd dimension for the moment !!\n";
                        // FIXME: mirror can only be done on the second dimension
                        if(tmp->numberOfDimensions()==2)
                        {
                            TH2F* histo = dynamic_cast<TH2F*>(tmp->getTemplate());
                            for (int binx=0;binx<histo->GetNbinsX(); binx++)
                            {
                                for (int biny=0;biny<histo->GetNbinsY()/2; biny++)
                                {
                                    double avr = histo->GetBinContent(binx+1,biny+1) - histo->GetBinContent(binx+1,histo->GetNbinsY()-biny);
                                    histo->SetBinContent(binx+1, biny+1, avr/2.);
                                    histo->SetBinContent(binx+1, histo->GetNbinsY()-biny, -avr/2.);
                                } 
                            }
                        }
                        else if(tmp->numberOfDimensions()==3)
                        {
                            TH3F* histo = dynamic_cast<TH3F*>(tmp->getTemplate());
                            for (int binx=0;binx<histo->GetNbinsX(); binx++)
                            {
                                for (int binz=0;binz<histo->GetNbinsZ(); binz++)
                                {
                                    for (int biny=0;biny<histo->GetNbinsY()/2; biny++)
                                    {
                                        double avr = histo->GetBinContent(binx+1,biny+1,binz+1) - histo->GetBinContent(binx+1,histo->GetNbinsY()-biny,binz+1);
                                        histo->SetBinContent(binx+1, biny+1, binz+1, avr/2.);
                                        histo->SetBinContent(binx+1, histo->GetNbinsY()-biny, binz+1, -avr/2.);
                                    } 
                                }
                            }
                        }
                        break;
                    }
                case Template::PostProcessing::FLOOR:
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
                        break;
                    }
                default:
                    break;
            }
            double sumOfWeightsAfter = tmp->getTemplate()->GetSumOfWeights();
            cout<<"[INFO]   Sum of weights after/before = "<<sumOfWeightsAfter<<" / "<<sumOfweightsBefore<<" = "<<sumOfWeightsAfter/sumOfweightsBefore<<"\n";
            sumOfweightsBefore = sumOfWeightsAfter;
        }
        // normalize
        //double sumOfWeights = tmp->getTemplate()->GetSumOfWeights();
        //tmp->getTemplate()->Scale(1./sumOfWeights);
    }
}


