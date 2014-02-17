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

#include <TH1.h>
#include "TCanvas.h"

#include <boost/any.hpp>

#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <map>
#include <stdexcept>



class PostProcessing
{
    public:
        enum Type
        {
            SMOOTH = 0,
            MIRROR = 1,
            FLOOR = 2,
            RESCALE = 3,
            REWEIGHT = 4
        };
        PostProcessing(Type type);
        ~PostProcessing();

        Type type() {return m_type;}
        template <class T> void addParameter(const std::string& par, const T& value)
        {
            m_parameters[par] = value;
        }

        template <class T> T getParameter(const std::string& par) const
        {
            std::map<std::string, boost::any>::const_iterator it = m_parameters.find(par);
            if(it==m_parameters.end())
            {
                std::stringstream error;
                error <<"PostProcessing::getParameter(): Trying to retrieve unknown parameter '"<<par<<"'\n";
                throw std::runtime_error(error.str());
            }
            boost::any value = it->second;
            if(value.type()!=typeid(T))
            {
                std::stringstream error;
                error <<"PostProcessing::getParameter(): Trying to retrieve parameter '"<<par<<"' with the wrong type\n";
                throw std::runtime_error(error.str());
            }
            return boost::any_cast<T>(value);
        }

    private:
        Type m_type;
        std::map<std::string, boost::any> m_parameters;
};

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
        unsigned int numberOfDimensions() const {return m_variables.size();} ;
        const std::string& getTreeName() const {return m_treeName;}
        double getRescaling() {return m_scaleFactor;}
        BinningType getBinningType() const {return m_binningType;}
        unsigned int getEntriesPerBin() const {return m_entriesPerBin;}
        TH1* getTemplate() const {return m_template;}
        TH1D* getRaw1DTemplate(unsigned int axis=0) const {return m_raw1DTemplates[axis];}
        const std::vector<TH1D*>& getRaw1DTemplates() const {return m_raw1DTemplates;}
        TH1D* getProjected1DTemplate(unsigned int axis=0);
        const std::vector< std::pair<double,double> >& getMinMax() const {return m_minmax;} 
        const std::vector<TH1*>& getWidths() const {return m_widths;}
        TH1* getWidth(int axis=0) const {return m_widths[axis];}
        std::vector< std::pair<std::string, std::string> >::const_iterator inputFileAndTreeBegin() const {return m_inputFileAndTreeNames.begin();}
        std::vector< std::pair<std::string, std::string> >::const_iterator inputFileAndTreeEnd() const {return m_inputFileAndTreeNames.end();}
        std::vector<std::pair<std::string,double> >::const_iterator inputTemplatesBegin() const {return m_inputTemplates.begin();}
        std::vector<std::pair<std::string,double> >::const_iterator inputTemplatesEnd() const {return m_inputTemplates.end();}
        std::vector<std::vector<double> >::const_iterator entriesBegin() const {return m_entries.begin();}
        std::vector<std::vector<double> >::const_iterator entriesEnd() const {return m_entries.end();}
        std::vector<double>::const_iterator weightsBegin() const {return m_weights.begin();}
        std::vector<double>::const_iterator weightsEnd() const {return m_weights.end();}
        const std::vector< std::vector<double> >& entries() const {return m_entries;}
        const std::vector<double>& weights() const {return m_weights;}
        std::vector<PostProcessing>::iterator postProcessingBegin() {return m_postProcessings.begin();}
        std::vector<PostProcessing>::iterator postProcessingEnd() {return m_postProcessings.end();}
        std::vector<TCanvas*>::iterator controlPlotsBegin() {return m_controlPlots.begin();}
        std::vector<TCanvas*>::iterator controlPlotsEnd() {return m_controlPlots.end();}

        // setters
        void setName(const std::string& name) {m_name = name;}
        void setOrigin(Origin origin) {m_origin = origin;}
        void setSelection(const std::string& selection) {m_selection = selection;}
        void setAssertion(const std::string& assertion) {m_assertion = assertion;}
        void setWeight(const std::string& weight) {m_weight = weight;}
        void setVariables(const std::vector<std::string>& vars);
        void addInputFileAndTree(const std::string& fileName, const std::string& treeName) {m_inputFileAndTreeNames.push_back( std::make_pair(fileName,treeName) );}
        void addInputTemplate(const std::string& name, double factor) {m_inputTemplates.push_back(make_pair(name,factor));}
        void setTreeName(const std::string& name);
        void setBinningType(BinningType type) {m_binningType = type;}
        void setEntriesPerBin(unsigned int entriesPerBin) {m_entriesPerBin = entriesPerBin;}
        void addPostProcessing(PostProcessing postProcess) {m_postProcessings.push_back(postProcess);}
        void createTemplate(const std::vector<unsigned int>& nbins, const std::vector< std::pair<double,double> >& minmax);
        void setTemplate(const TH1* histo);
        void setRaw1DTemplates(const std::vector<TH1D*>& histo);
        void setWidths(const std::vector<TH1*>& width);
        void setRescaling(double scaleFactor) {m_scaleFactor = scaleFactor;}
        void store(const std::vector<double>& vs, double w);
        void reweight1D(unsigned int axis, unsigned int bin, double weight);
        // control plot methods
        void makeProjectionControlPlot(const std::string& tag);
        void addControlPlot(TCanvas* plot) {m_controlPlots.push_back(plot);}


    private:
        std::string m_name;
        Origin m_origin;
        std::vector< std::pair<std::string,std::string> > m_inputFileAndTreeNames;
        std::vector<std::pair<std::string, double> > m_inputTemplates;
        std::string m_treeName;
        std::string m_selection;
        std::string m_assertion;
        std::string m_weight;
        std::vector<std::string> m_variables;
        BinningType m_binningType;
        TH1* m_template;
        std::vector<TH1D*> m_raw1DTemplates;
        std::vector< std::pair<double,double> > m_minmax;
        std::vector<TH1*> m_widths;
        unsigned int m_entriesPerBin;
        std::vector<PostProcessing> m_postProcessings;
        double m_scaleFactor;
        std::vector< std::vector<double> > m_entries;
        std::vector< double > m_weights;

        std::vector<TCanvas*> m_controlPlots;

};


#endif
