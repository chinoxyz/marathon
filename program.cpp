#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <set>
#include <queue>
#include <stack>
#include <string>
#include <sstream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cassert>

using namespace std;

#define ri(X) scanf("%d", &(X))
#define pi(X) printf("%d", (X))
#define mp(X,Y) make_pair(X,Y)
#define pb(X) push_back(X)
#define lint long long
#define pii pair<int,int>
#define inf 1e9
#define linf 1e18
#define X first
#define Y second
#define all(X) (X).begin(),(X).end()
#define uni(X) X.erase(unique(X.begin(), X.end()), X.end());


class PT{
  public:
  int x, y, i;
  PT(int xx = 0, int yy = 0, int ii = 0){
    x = xx;
    y = yy;
    i = ii;
  }
};

string itos(int x){
  stringstream ss;
  ss << x;
  return ss.str();
}

class SmallPolygons{
  vector<PT> ve;
  vector<string> res;
  public:
  vector <string> choosePolygons(vector <int> points, int N){
    int np = points.size();
    
    for(int i = 0; i < np; i+=2){
      ve.pb(PT(points[i],points[i+1],i/2));
    }
    np/=2;
    string rp;
    //comment
    for(int i = 0; i < np; i++){
      rp+=itos(i)+ " ";
    }
    res.pb(rp);
    return res;
  }

};
int main(){
  int Np; 
  vector<int> points;
  int N;
  int v;
  vector<string> ret;
  cin >> Np;
  for (int i=0; i < Np; i++){
    cin >> v;
    points.pb(v);
  }
  cin >> N;
  
  ret = SmallPolygons().choosePolygons(points, N);
  printf("%d\n", ret.size());
  for (int i=0; i < ret.size(); i++){
    printf("%s\n", ret[i].c_str());
    fprintf(stderr,"%s\n", ret[i].c_str());
  }
  fflush(stdout);

  return 0;
}
