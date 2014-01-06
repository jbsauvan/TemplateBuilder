/**
 *  @file  main.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    12/09/2013
 *
 *  @internal
 *     Created :  12/09/2013
 * Last update :  12/09/2013 11:10:29 PM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */



#include <string>
#include <iostream>
#include <stdexcept>


#include "TemplateManager.h"

int main(int argc, char** argv)
{
    if(argc!=2)
    {
        std::cerr<<"Usage: buildtemplate.exe parFile.json\n";
        return EXIT_FAILURE;
    }

    std::string parFile(argv[1]);

    TemplateManager manager;
    try
    {
        manager.initialize(parFile);
        manager.loop();
    }catch(std::exception& e)
    {
        std::cerr<<"[ERROR] "<<e.what()<<"\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
