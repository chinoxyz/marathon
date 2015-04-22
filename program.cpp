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



string itos(int x){
  stringstream ss;
  ss << x;
  return ss.str();
}

string vtos(vector<int> vi){
  string rp;
  for(int i = 0; i < vi.size(); i++){
    rp+=itos(vi[i])+ " ";
  }
  return rp;
}
vector<string> vvtovs(vector<vector<int> > res){
  vector<string> vs;
  for(int i = 0; i < res.size(); i++){
    vs.pb(vtos(res[i]));
  }
  return vs;
  
}



const double pi = acos(-1);
class PT{
  public:
  int x, y, i;
  PT(int xx = 0, int yy = 0, int ii = 0){
    x = xx;
    y = yy;
    i = ii;
  }
  bool operator <(const PT &p)  const { return x < p.x || (x == p.x && y < p.y); }
  bool operator ==(const PT &p)  const { return x == p.x && y == p.y; }
  PT operator + (const PT &p)  const { return PT(x+p.x, y+p.y); }
  PT operator - (const PT &p)  const { return PT(x-p.x, y-p.y); }
  PT operator * (double c)     const { return PT(x*c,   y*c  ); }
  PT operator / (double c)     const { return PT(x/c,   y/c  ); }
};
int dot(PT p, PT q)     { return p.x*q.x+p.y*q.y; }
int dist2(PT p, PT q)   { return dot(p-q,p-q); }
double dist(PT p, PT q)   { return sqrt(dist2(p,q)); }
int cross(PT p, PT q)   { return p.x*q.y-p.y*q.x; }
double angle(PT o, PT p, PT q){
  return ((cross(p-o,q-o)<0?-1:1)*acos(dot(p-o,q-o)/(dist(o,p)*dist(o,q))));
}

double angle2(PT o, PT p, PT q){
  return asin(cross(p-o,q-o)/(dist(o,p)*dist(o,q)));
}

ostream &operator<<(ostream &os, const PT &p) {
  os << "(" << p.x << "," << p.y << ")"; 
}

vector<PT> ve;
vector<string> res;
int np;
int n;


class solution{
public:
  vector<vector<int> > ve;
};

vector<int> getpol1(vector<int> vi){
  assert(vi.size() >= 3);
  if(vi.size() == 3) return vi;
  vector<PT> vp;
  for(int i = 0;i < vi.size(); i++){
    vp.pb(ve[vi[i]]);
  }
  sort(all(vp));
  
  PT ic = vp[0];
  PT ir = ic+PT(-1000,0);
  
  
  vector<pair<pair<double,double>, int> > vs;
  //cerr << ic << endl;
  
  for(int i = 0; i < vi.size(); i++){
    if(ve[vi[i]] == ic){
      vs.pb(mp(mp(-1,-1),i));
      continue;
    } 
    double ang = angle(ic, ve[vi[i]], ir);
    double di = dist2(ic, ve[vi[i]]);
    vs.pb(mp(mp(ang,di),i));
  }
  sort(all(vs));
  for(int i = 0; i < vs.size(); i++){
    vi[i] = vs[i].Y;
  }
  /*
  for(int i = vs.size(); i >= 0; i--){
    cerr << ve[vs[i].Y] << " " << vs[i].X.X << " " << vs[i].X.Y << endl;
  }*/
  return vi;
}

vector<vector<int> > getsol1(){
  vector<vector<int> > res;
  vector<int> vi;
  for(int i = 0; i < np; i++){
    vi.pb(i);
  }
  res.pb(vi);
  for(int i = 0; i < res.size(); i++){
    res[i] = getpol1(res[i]);
  }
  return res;
}

vector<vector<int> > getsol2(){
return getsol1();
  if(n == 1) return getsol1();
  
  vector<vector<int> > res;
  vector<PT> vp = ve;
  sort(all(vp));
  
  vector<int> vb;
  for(int i = 0; i < n-1; i++){
    vb.pb(rand()%np);
  }vb.pb(np-1);
  
  sort(all(vb));
  uni(vb);
  reverse(all(vb));
  
  
  for(int i = 0; i < 2; i++){
    res.pb(vector<int>());
  }
  
  for(int i = 0; i < np; i++){
    res[(i<np/2?0:1)].pb(i);//vp[i].i
    /*vi.pb(i);//vp[i].i
    fprintf(stderr,"%d ",vp[i].i);
    if(!vb.empty() && vb.back() == i && vi.size() >= 3){
      res.pb(vi);
      vi.clear();
      vb.pop_back();
      fprintf(stderr,"\n");
    }*/
  }
  for(int i = 0; i < res.size(); i++){
    res[i] = getpol1(res[i]);
  }
  return res;
}


class SmallPolygons{
  
    
  
  
  
  
  
  public:
  vector <string> choosePolygons(vector <int> points, int N){
    np = points.size();
    n = N;
    for(int i = 0; i < np; i+=2){
      ve.pb(PT(points[i],points[i+1],i/2));
    }
    np/=2;
    
    
    
    
    vector<vector<int> > res = getsol2();
    
    return vvtovs(res);
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
    //fprintf(stderr,"VEC: %s\n", ret[i].c_str());
  }
  fflush(stdout);

  return 0;
}
