#include<iostream>
#include<set>
#include<vector>
#include<algorithm>

using namespace std;

const int maxn = 2000000;

int n, m;
vector<pair<pair<int, int>, int> > all;
set<pair<int, int>, greater<pair<int, int> > > S;
pair<int, int> A[maxn], B[maxn];
pair<int, int> ans[maxn];

int main(){
	cin >> n >> m;
	for(int i = 0; i < n; i++){
		cin >> A[i].first >> A[i].second;
		all.push_back(make_pair(make_pair(A[i].first, i), 0));
		all.push_back(make_pair(make_pair(A[i].second, i), 2));
	}
	for(int i = 0; i < m; i++){
		cin >> B[i].first >> B[i].second;
		all.push_back(make_pair(make_pair(B[i].first, i), 1));
	}
	sort(all.begin(), all.end());
	for(int i = 0; i < all.size(); i++){
		int x = all[i].first.first;
		int ind = all[i].first.second;
		int t = all[i].second;
		if(t == 0){
			S.insert(make_pair(A[ind].second, ind));
		}
		else if(t == 2){
			S.erase(make_pair(x, ind));
		}
		else{
			auto mx = S.begin()->first;
			auto interval = S.begin()->second;
			if(mx >= B[ind].second){
				ans[ind] = A[interval];
			}
			else{
				ans[ind] = make_pair(-1, -1);
			}
		}
	}

	for(int i = 0; i < m; i++){
		if (ans[i] != make_pair(-1,-1)){
cout  << B[i].first << "\t" << B[i].second << "\t" << ans[i].first << "\t" << ans[i].second  << endl; 	}
	}
	return 0;
}

