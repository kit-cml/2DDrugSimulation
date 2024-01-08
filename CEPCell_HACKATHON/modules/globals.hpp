#ifndef GLOBALS_HPP
#define GLOBALS_HPP

class mympi
{
 public:
    static int rank, size;
};

class glob_var
{
  public:
    static int A1656D_mode; // 1: Mutation, 2: Mexilatine, 3: Flecaindine, 4: Ranolazine
    static bool is_hrv;
};

#endif
