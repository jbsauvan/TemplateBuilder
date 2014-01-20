/**
 *  @file  TemplateManager.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    12/08/2013
 *
 *  @internal
 *     Created :  12/08/2013
 * Last update :  12/08/2013 02:45:24 PM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */

#include "TemplateManager.h"

#include <TTree.h>
#include <TFile.h>
#include <TEntryList.h>
#include <TTreeFormula.h>

#include <math.h> 
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>


using namespace std;

/*****************************************************************/
TemplateManager::TemplateManager():
    m_inputDirectory("./"),
    m_outputFileName("templates.root"),
    m_outputFile(NULL)
/*****************************************************************/
{
}


/*****************************************************************/
TemplateManager::~TemplateManager()
/*****************************************************************/
{
    if(m_outputFile)
    {
        m_outputFile->Write();
        m_outputFile->Close();
    }
}


/*****************************************************************/
bool TemplateManager::initialize(const std::string& parameterFile)
/*****************************************************************/
{
    cout<<"[INFO] Initializing\n";
    m_reader.read(parameterFile);

    m_inputDirectory = m_reader.inputDirectory();
    struct stat buf;
    stat(m_inputDirectory.c_str(), &buf);
    if(!S_ISDIR(buf.st_mode))
    {
        stringstream error;
        error <<"TemplateManager::initialize(): Input directory '"<<m_inputDirectory<<"' doesn't exist\n";
        throw runtime_error(error.str());
    }
    cout<<"[INFO]   Input directory: "<<m_inputDirectory<<"\n";


    m_outputFileName = m_reader.outputFileName();
    m_outputFile = TFile::Open(m_outputFileName.c_str(), "RECREATE");
    if(!m_outputFile)
    {
        stringstream error;
        error <<"Cannot open output file for writing\n";
        throw runtime_error(error.str());
    }
    cout<<"[INFO]   Output file: "<<m_outputFileName<<"\n";

    unsigned int nTemplates = 0;
    vector<Template*>::iterator tmpIt = m_reader.templateBegin();
    vector<Template*>::iterator tmpItE = m_reader.templateEnd();
    for(;tmpIt!=tmpItE;++tmpIt)
    {
        m_templates.addTemplate(*tmpIt);
        nTemplates++;
    }
    cout<<"[INFO]   "<<nTemplates<<" templates will be produced\n";


    return true;
}


/*****************************************************************/
void TemplateManager::loop()
/*****************************************************************/
{
    map<string,Template*>::iterator tmpIt = m_templates.templateBegin();
    map<string,Template*>::iterator tmpItE = m_templates.templateEnd();
    for(;tmpIt!=tmpItE;++tmpIt)
    {
        const string& tmpName = tmpIt->first;
        Template* tmp = tmpIt->second;
        if(tmp->getOrigin()!=Template::Origin::FILES) continue;

        cout<<"[INFO] Filling template '"<<tmpName<<"' from trees\n";
        if(tmp->getWeight()!="")
        {
            cout<<"[INFO]   Weights '"<<tmp->getWeight()<<"' will be used\n";
        }
        vector<string>::const_iterator fIt = tmp->inputFileBegin();
        vector<string>::const_iterator fItE = tmp->inputFileEnd();
        int nFiles = (int)(fItE-fIt);
        int i = 1;
        for(;fIt!=fItE;++fIt)
        {
            const string& fileName = *fIt;
            stringstream fullName;
            fullName << m_inputDirectory<< "/" << fileName;
            std::cout<<"[INFO]   Opening file "<<i<<"/"<<nFiles<<"\n";
            TFile* inputFile = TFile::Open(fullName.str().c_str());
            if(!inputFile)
            {
                stringstream error;
                error << "TemplateManager::loop(): Cannot open file '"<<fileName<<"'\n";
                throw runtime_error(error.str());
            }
            TTree* tree = dynamic_cast<TTree*>(inputFile->Get(tmp->getTreeName().c_str()));
            tree->Draw(">>elist", tmp->getSelection().c_str(), "entrylist");
            TEntryList* entryList = dynamic_cast<TEntryList*>(gDirectory->Get("elist"));
            Long64_t nEntries = entryList->GetN();
            tree->SetEntryList(entryList);

            vector<TTreeFormula*> varForms;
            for(unsigned int v=0;v<tmp->numberOfDimensions();v++)
            {
                stringstream varName;
                varName << "var" << v;
                varForms.push_back( new TTreeFormula(varName.str().c_str(), tmp->getVariable(v).c_str(), tree) );
            }
            TTreeFormula* weightForm = NULL;
            if(tmp->getWeight()!="")
            {
                weightForm = new TTreeFormula("weight", tmp->getWeight().c_str(), tree);
            }
            TTreeFormula assertForm("assert", tmp->getAssertion().c_str(), tree);

            // first compute sum of weights for global normalization of the sample
            double sumOfWeights = (double)nEntries;
            if(weightForm)
            {
                sumOfWeights = 0.;
                for (Long64_t entry=0;entry<nEntries;entry++)
                {
                    tree->GetEntry(entry);
                    double weight = weightForm->EvalInstance();
                    if(!std::isfinite(weight))
                    {
                        std::cerr<<"[WARNING]   Inf or NaN weight\n";
                    }
                    sumOfWeights += weight;
                }
                if(sumOfWeights==0)
                {
                    std::cerr<<"[WARNING]   Sum of weights = 0\n";
                }
                //cout<<"Computing sum of weights: "<<nEntries<<"/"<<sumOfWeights<<"\n";
            }

            for (Long64_t entry=0;entry<nEntries;entry++)
            {
                //if(entry%10000==0)
                //{
                //    cout<<"[INFO] Getting entry "<<entry<<"/"<<nEntries<<"\n";
                //}

                Long64_t entryNumber = tree->GetEntryNumber(entry); 
                tree->GetEntry(entryNumber);
                if(!assertForm.EvalInstance())
                {
                    stringstream error;
                    error << "TemplateManager::loop(): assertion '"<<tmp->getAssertion()<<"' failed";
                    throw runtime_error(error.str());
                }
                vector<double> point;
                for(unsigned int v=0;v<tmp->numberOfDimensions();v++)
                {
                    double varValue = varForms[v]->EvalInstance();
                    if(!std::isfinite(varValue))
                    {
                        std::cerr<<"[WARNING]   Inf or NaN variable\n";
                    }
                    point.push_back(varValue);
                }
                double weight = 1.;
                if(weightForm)
                {
                    weight = weightForm->EvalInstance();
                }
                tmp->store(point, weight*(double)nEntries/sumOfWeights);

                //execute();
            }
            if(weightForm)
            {
                weightForm->Delete();
            }
            for(unsigned int v=0;v<tmp->numberOfDimensions();v++)
            {
                varForms[v]->Delete();
            }
            inputFile->Close();
            i++;
        }
    }
    m_templates.fillTemplates();
    m_templates.postProcessing();
    m_templates.buildTemplatesFromTemplates();
    save();
}


/*****************************************************************/
void TemplateManager::save()
/*****************************************************************/
{
    m_outputFile->cd();
    map<string,Template*>::iterator tmpIt = m_templates.templateBegin();
    map<string,Template*>::iterator tmpItE = m_templates.templateEnd();
    for(;tmpIt!=tmpItE;++tmpIt)
    {
        const string& tmpName = tmpIt->first;
        Template* tmp = tmpIt->second;
        tmp->getTemplate()->SetName(tmpName.c_str());
        tmp->getTemplate()->Write();

        // TMP: fill kernel widths
        //tmp->getWidth(0)->Write();
        //tmp->getWidth(1)->Write();
    }
}

