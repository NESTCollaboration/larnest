/**
 * @file Logger.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-06-27
 */
#include "Logger.h"

namespace larnest 
{
    Logger* Logger::sInstance{nullptr};
    std::mutex Logger::sMutex;

    Logger *Logger::GetInstance(const std::string& name)
    {
        std::lock_guard<std::mutex> lock(sMutex);
        if (sInstance == nullptr)
        {
            sInstance = new Logger(name);
        }
        sInstance->SetName(name);
        return sInstance;
    }

    Logger::Logger(const std::string name)
    : sName(name)
    {
        std::time_t now = std::time(0);
        std::tm *ltm = std::localtime(&now);
        sYear = std::to_string(1900 + ltm->tm_year);
        sMonth = std::to_string(1 + ltm->tm_mon);
        sDay = std::to_string(ltm->tm_mday);
        sDate = sMonth + "-" + sDay + "-" + sYear;
        mkdir(".logs", 0777);
        sOutputFileName = ".logs/larnest_" + sDate + ".log";
    }
}