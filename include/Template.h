/**
 *  @file  Template.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    12/08/2013
 *
 *  @internal
 *     Created :  12/08/2013
 * Last update :  12/08/2013 02:54:03 PM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */



#ifndef TEMPLATE_H
#define TEMPLATE_H

#include <TH2F.h>

#include <string>
#include <vector>
#include <assert.h>
#include <iostream>


class Template
{
    public:
        enum BinningType
        {
            FIXED = 0,
            ADAPTIVE = 1
        };
        enum Origin
        {
            FILES = 0,
            TEMPLATES = 1
        };
        enum PostProcessing
        {
            SMOOTH_K5B = 0,
            SMOOTH_ADAPTIVE = 1,
            MIRROR = 2,
            MIRROR_INV = 3,
            FLOOR = 4
        };


        Template();
        Template(const Template& tmp);
        ~Template();

        // getters
        const std::string& getName()const {return m_name;}
        Origin getOrigin() const {return m_origin;}
        const std::string& getSelection() const {return m_selection;}
        const std::string& getAssertion() const {return m_assertion;}
        const std::string& getWeight() const {return m_weight;}
        const std::string& getVariable(unsigned int index) const;
        const std::string& getTreeName() const {return m_treeName;}
        double getRescaling() {return m_scaleFactor;}
        BinningType getBinningType() const {return m_binningType;}
        unsigned int getEntriesPerBin() const {return m_entriesPerBin;}
        TH2F* getTemplate() const {return m_template;}
        TH2F* getWidthX() const {return m_widthx;}
        TH2F* getWidthY() const {return m_widthy;}
        std::vector<std::string>::const_iterator inputFileBegin() const {return m_inputFileNames.begin();}
        std::vector<std::string>::const_iterator inputFileEnd() const {return m_inputFileNames.end();}
        std::vector<std::pair<std::string,double> >::const_iterator inputTemplatesBegin() const {return m_inputTemplates.begin();}
        std::vector<std::pair<std::string,double> >::const_iterator inputTemplatesEnd() const {return m_inputTemplates.end();}
        std::vector<std::vector<double> >::const_iterator entriesBegin() const {return m_entries.begin();}
        std::vector<std::vector<double> >::const_iterator entriesEnd() const {return m_entries.end();}
        const std::vector<std::vector<double> >& entries() const {return m_entries;}
        std::vector<PostProcessing>::iterator postProcessingBegin() {return m_postProcessings.begin();}
        std::vector<PostProcessing>::iterator postProcessingEnd() {return m_postProcessings.end();}

        // setters
        void setName(const std::string& name) {m_name = name;}
        void setOrigin(Origin origin) {m_origin = origin;}
        void setSelection(const std::string& selection) {m_selection = selection;}
        void setAssertion(const std::string& assertion) {m_assertion = assertion;}
        void setWeight(const std::string& weight) {m_weight = weight;}
        void setVariables(const std::string& var1, const std::string& var2);
        void addInputFile(const std::string& name) {m_inputFileNames.push_back(name);}
        void addInputTemplate(const std::string& name, double factor) {m_inputTemplates.push_back(make_pair(name,factor));}
        void setTreeName(const std::string& name) {m_treeName = name;}
        void setBinningType(BinningType type) {m_binningType = type;}
        void setEntriesPerBin(unsigned int entriesPerBin) {m_entriesPerBin = entriesPerBin;}
        void addPostProcessing(PostProcessing postProcess) {m_postProcessings.push_back(postProcess);}
        void createTemplate(unsigned int nbin1, double min1, double max1, unsigned int nbin2, double min2, double max2);
        void setTemplate(const TH2F* histo);
        void setWidths(const TH2F* widthx, const TH2F* widthy);
        void setRescaling(double scaleFactor) {m_scaleFactor = scaleFactor;}
        void store(double v1, double v2, double w);


    private:
        std::string m_name;
        Origin m_origin;
        std::vector<std::string> m_inputFileNames;
        std::vector<std::pair<std::string, double> > m_inputTemplates;
        std::string m_treeName;
        std::string m_selection;
        std::string m_assertion;
        std::string m_weight;
        std::vector<std::string> m_variables;
        BinningType m_binningType;
        TH2F* m_template;
        TH2F* m_widthx;
        TH2F* m_widthy;
        unsigned int m_entriesPerBin;
        std::vector<PostProcessing> m_postProcessings;
        double m_scaleFactor;
        std::vector< std::vector<double> > m_entries;

};


#endif
