
#include "TemplateBuilder.h"
#include "Bin2DTree.h"
#include "GaussKernelSmoother.h"

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
        error << "TemplateBuilder::getTemplate(): template " <<name<<" doesnt' exist";
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
            cout<< "[INFO] Building template '"<<tmp->getName()<<"' with standard binning\n";
            vector<vector<double> >::const_iterator it  = tmp->entriesBegin();
            vector<vector<double> >::const_iterator itE = tmp->entriesEnd();
            for(;it!=itE;++it)
            {
                vector<double> entry = *it;
                tmpIt->second->getTemplate()->Fill(entry[0],entry[1],entry[2]);
            }
        }
        else if(tmp->getBinningType()==Template::BinningType::ADAPTIVE)
        {
            cout<< "[INFO] Deriving adaptive binning for template '"<<tmp->getName()<<"'... This may take some time\n";
            double min1 = tmp->getTemplate()->GetXaxis()->GetBinLowEdge(1);
            double max1 = tmp->getTemplate()->GetXaxis()->GetBinUpEdge(tmp->getTemplate()->GetNbinsX());
            double min2 = tmp->getTemplate()->GetYaxis()->GetBinLowEdge(1);
            double max2 = tmp->getTemplate()->GetYaxis()->GetBinUpEdge(tmp->getTemplate()->GetNbinsY());
            Bin2DTree bintree(min1, max1, min2, max2, tmp->entries());
            bintree.setMinLeafEntries(tmp->getEntriesPerBin());
            TH2F* gridConstraint = (TH2F*)tmp->getTemplate()->Clone("gridConstraint");
            bintree.setGridConstraint(gridConstraint);
            bintree.build();
            cout<<"[INFO]   Number of bins = "<<bintree.getNLeaves()<<"\n";
            TH2F* histo = dynamic_cast<TH2F*>(bintree.fillHistogram());
            tmp->setTemplate(histo);
            pair<TH2F*,TH2F*> widths = bintree.fillWidths();
            tmp->setWidths(widths.first, widths.second);
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
                        cout<<"[INFO] Smoothing template '"<<tmp->getName()<<"' with k5b kernel\n";
                        tmp->getTemplate()->Smooth(1, "k5b");
                        break;
                    }
                case Template::PostProcessing::SMOOTH_ADAPTIVE:
                    {
                        cout<<"[INFO] Smoothing template '"<<tmp->getName()<<"' with variable Gaussian kernel... This may take some time depending on the number of bins\n";
                        GaussKernelSmoother smoother;
                        smoother.setWidths(tmp->getWidthX(), tmp->getWidthY());
                        TH2F* histoSmooth = smoother.smooth(tmp->getTemplate());
                        tmp->setTemplate(histoSmooth);
                        break;
                    }
                case Template::PostProcessing::MIRROR:
                    {
                        cout<<"[INFO] Mirroring template '"<<tmp->getName()<<"'\n";
                        TH2F* histo = tmp->getTemplate();
                        for (int binx=0;binx<histo->GetNbinsX(); binx++)
                        {
                            for (int biny=0;biny<histo->GetNbinsY()/2; biny++)
                            {
                                double avr = histo->GetBinContent(binx+1,biny+1) + histo->GetBinContent(binx+1,histo->GetNbinsY()-biny);
                                histo->SetBinContent(binx+1, biny+1, avr/2.);
                                histo->SetBinContent(binx+1, histo->GetNbinsY()-biny, avr/2.);
                            } 
                        }
                        break;
                    }
                case Template::PostProcessing::MIRROR_INV:
                    {
                        cout<<"[INFO] Anti-mirroring template '"<<tmp->getName()<<"'\n";
                        TH2F* histo = tmp->getTemplate();
                        for (int binx=0;binx<histo->GetNbinsX(); binx++)
                        {
                            for (int biny=0;biny<histo->GetNbinsY()/2; biny++)
                            {
                                double avr = histo->GetBinContent(binx+1,biny+1) - histo->GetBinContent(binx+1,histo->GetNbinsY()-biny);
                                histo->SetBinContent(binx+1, biny+1, avr/2.);
                                histo->SetBinContent(binx+1, histo->GetNbinsY()-biny, -avr/2.);
                            } 
                        }
                        break;
                    }
                case Template::PostProcessing::FLOOR:
                    {
                        cout<<"[INFO] Flooring template '"<<tmp->getName()<<"'\n";
                        TH2F* histo = tmp->getTemplate();
                        double floorN = ((histo->Integral())/(histo->GetNbinsX()*histo->GetNbinsY()))*(0.001/100.);
                        for(int binx = 1; binx <= histo->GetNbinsX(); binx++)
                        {
                            for(int biny = 1; biny <= histo->GetNbinsY(); biny++)
                            {
                                double orig = histo->GetBinContent(binx,biny);
                                histo->SetBinContent(binx,biny,(orig+floorN));
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
                tmp->setWidths(inTmp->getWidthX(), inTmp->getWidthY());
            }
            else
            {
                // TODO: take the averaged width for each bin instead of the width of the first template
                cout<<"[INFO] + ("<<factor<<") x "<<inTmp->getName()<<"\n";
                int nbinsx1 = tmp->getTemplate()->GetNbinsX();
                int nbinsy1 = tmp->getTemplate()->GetNbinsY();
                double minx1 = tmp->getTemplate()->GetXaxis()->GetBinLowEdge(0);
                double maxx1 = tmp->getTemplate()->GetXaxis()->GetBinUpEdge(nbinsx1);
                double miny1 = tmp->getTemplate()->GetYaxis()->GetBinLowEdge(0);
                double maxy1 = tmp->getTemplate()->GetYaxis()->GetBinUpEdge(nbinsy1);
                int nbinsx2 = inTmp->getTemplate()->GetNbinsX();
                int nbinsy2 = inTmp->getTemplate()->GetNbinsY();
                double minx2 = inTmp->getTemplate()->GetXaxis()->GetBinLowEdge(0);
                double maxx2 = inTmp->getTemplate()->GetXaxis()->GetBinUpEdge(nbinsx2);
                double miny2 = inTmp->getTemplate()->GetYaxis()->GetBinLowEdge(0);
                double maxy2 = inTmp->getTemplate()->GetYaxis()->GetBinUpEdge(nbinsy2);
                if(nbinsx1!=nbinsx2 || nbinsy1!=nbinsy2 ||
                        minx1!=minx2 || maxx1!=maxx2 ||
                        miny1!=miny2 || maxy1!=maxy2
                        )
                {
                    stringstream error;
                    error << "TemplateBuilder::buildTemplatesFromTemplates(): Trying to add templates with different binning ("<<tmp->getName()<<" & "<<inTmp->getName()<<")\n";
                    throw runtime_error(error.str());
                }
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

        vector<Template::PostProcessing>::iterator it = tmp->postProcessingBegin();
        vector<Template::PostProcessing>::iterator itE = tmp->postProcessingEnd();
        for(;it!=itE;++it)
        {
            switch(*it)
            {
                case Template::PostProcessing::SMOOTH_K5B:
                    {
                        cout<<"[INFO] Smoothing template '"<<tmp->getName()<<"' with k5b kernel\n";
                        tmp->getTemplate()->Smooth(1, "k5b");
                        break;
                    }
                case Template::PostProcessing::SMOOTH_ADAPTIVE:
                    {
                        cout<<"[INFO] Smoothing template '"<<tmp->getName()<<"' with variable Gaussian kernel... This may take some time depending on the number of bins\n";
                        GaussKernelSmoother smoother;
                        smoother.setWidths(tmp->getWidthX(), tmp->getWidthY());
                        TH2F* histoSmooth = smoother.smooth(tmp->getTemplate());
                        tmp->setTemplate(histoSmooth);
                        break;
                    }
                case Template::PostProcessing::MIRROR:
                    {
                        cout<<"[INFO] Mirroring template '"<<tmp->getName()<<"'\n";
                        TH2F* histo = tmp->getTemplate();
                        for (int binx=0;binx<histo->GetNbinsX(); binx++)
                        {
                            for (int biny=0;biny<histo->GetNbinsY()/2; biny++)
                            {
                                double avr = histo->GetBinContent(binx+1,biny+1) + histo->GetBinContent(binx+1,histo->GetNbinsY()-biny);
                                histo->SetBinContent(binx+1, biny+1, avr/2.);
                                histo->SetBinContent(binx+1, histo->GetNbinsY()-biny, avr/2.);
                            } 
                        }
                        break;
                    }
                case Template::PostProcessing::MIRROR_INV:
                    {
                        cout<<"[INFO] Anti-mirroring template '"<<tmp->getName()<<"'\n";
                        TH2F* histo = tmp->getTemplate();
                        for (int binx=0;binx<histo->GetNbinsX(); binx++)
                        {
                            for (int biny=0;biny<histo->GetNbinsY()/2; biny++)
                            {
                                double avr = histo->GetBinContent(binx+1,biny+1) - histo->GetBinContent(binx+1,histo->GetNbinsY()-biny);
                                histo->SetBinContent(binx+1, biny+1, avr/2.);
                                histo->SetBinContent(binx+1, histo->GetNbinsY()-biny, -avr/2.);
                            } 
                        }
                        break;
                    }
                case Template::PostProcessing::FLOOR:
                    {
                        cout<<"[INFO] Flooring template '"<<tmp->getName()<<"'\n";
                        TH2F* histo = tmp->getTemplate();
                        double floorN = ((histo->Integral())/(histo->GetNbinsX()*histo->GetNbinsY()))*(0.001/100.);
                        for(int binx = 1; binx <= histo->GetNbinsX(); binx++)
                        {
                            for(int biny = 1; biny <= histo->GetNbinsY(); biny++)
                            {
                                double orig = histo->GetBinContent(binx,biny);
                                histo->SetBinContent(binx,biny,(orig+floorN));
                            }
                        }
                        break;
                    }
                default:
                    break;
            }
        }
        // normalize
        //double sumOfWeights = tmp->getTemplate()->GetSumOfWeights();
        //tmp->getTemplate()->Scale(1./sumOfWeights);
    }
}


