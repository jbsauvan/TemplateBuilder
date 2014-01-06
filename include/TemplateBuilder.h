/**
 *  @file  TemplateBuilder.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    12/08/2013
 *
 *  @internal
 *     Created :  12/08/2013
 * Last update :  12/08/2013 03:07:51 PM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */



#ifndef TEMPLATEBUILDER_H
#define TEMPLATEBUILDER_H

#include "Template.h"

#include <map>
#include <string>

class TemplateBuilder
{
    public:
        TemplateBuilder(){};
        ~TemplateBuilder();

        void addTemplate(const std::string& name);
        void addTemplate(Template* tmp);
        Template* getTemplate(const std::string& name);
        std::map<std::string, Template*>::const_iterator templateBegin() const {return m_templates.begin();}
        std::map<std::string, Template*>::const_iterator templateEnd() const {return m_templates.end();}
        std::map<std::string, Template*>::iterator templateBegin() {return m_templates.begin();}
        std::map<std::string, Template*>::iterator templateEnd() {return m_templates.end();}

        void fillTemplates();
        void postProcessing();
        void buildTemplatesFromTemplates();


    private:
        std::map<std::string, Template*> m_templates;
};


#endif
