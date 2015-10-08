

namespace std{


class Orientational_Correlation: public MULTIBODY_ANALYSIS
{
    float * correlation;
    float * weighting;
    float overall_correlation
    int n_times;
    Coordinate correlated_vector;

public:
    Orientational_Correlation();
    Orientational_Correlation(const Orientational_Correlation &);
    Orientational_Correlation operator = (const Orientational_Correlation &);

    Orientational_Correlation(System*);
    Orientaitonal_Correlation(System*, Coordinate);

    void analyze(Multibody_List * mblist);
    void list_displacementkernel(int,int,int){};
    void listkernel(Multibody *, int, int, int);
    void postprocess_list();

    void write(string) const;

};

}
