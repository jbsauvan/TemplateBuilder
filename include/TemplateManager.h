/**
 *  @file  TemplateManager.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    12/08/2013
 *
 *  @internal
 *     Created :  12/08/2013
 * Last update :  12/08/2013 02:43:44 PM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */



#ifndef TEMPLATEMANAGER_H
#define TEMPLATEMANAGER_H

#include "TemplateBuilder.h"
#include "TemplateParameters.h"

#include <string>

class TTree;
class TEntryList;
class TFile;

class TemplateManager
{
    public:
        TemplateManager();
        ~TemplateManager();

        bool initialize(const std::string& parameterFile);
        void loop();
        void fillTemplate();
        void save();


    protected:
        std::string m_inputDirectory;
        std::string m_outputFileName;
        TTree* m_inputTree;
        TEntryList* m_entryList;
        //Long64_t m_nEntries;
        TFile* m_outputFile;

        TemplateBuilder m_templates;
        TemplateParameters m_reader;

};

#endif
