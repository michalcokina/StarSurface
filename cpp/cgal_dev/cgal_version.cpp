#include <CGAL/config.h>
#include <iostream>
#include <iostream>
#include <sstream>
#include <string>

std::string IntToString ( int number ) {
  std::ostringstream oss;
  // Works just like cout
  oss<<number;
  // Return the underlying string
  return oss.str();
}

int StringToInt (std::string s) {
    int i = 0;
    std::istringstream ss;
    ss.str(s);
    ss >> i;
    return i;
}

float StringToFloat (std::string s) {
    float i;
    std::istringstream ss;
    ss.str(s);
    ss >> i;
    return i;
}

int main() {
    // cgal version
    int v = CGAL_VERSION_NR;
    std::string cgal_version = IntToString(v);
    std::string main_v = cgal_version.substr(1, 2);
    std::string minor_v = cgal_version.substr(3, 2);
    std::string version_s = IntToString(StringToInt(main_v)) + "." + IntToString(StringToInt(minor_v));
    float version = StringToFloat(version_s);
    
    return 0;
}
