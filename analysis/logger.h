#ifndef scribe_H
#define scribe_H

namespace scribe {
  class logType{};
  class DebugLog : logType{};
  class NormalLog : logType{};
  
  enum class LogLevel{
    low,
      medium,
      high
      };

  class Scribe{
    private:
    LogLevel debugLog(DebugLog){return }
    template<typename Object, typename Method, LogLevel l = LogLevel::low> void log(Object o, Method m){
      if()


}

  };
}
#endif
