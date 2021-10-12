#include <iostream>
#include <string>
#include <unordered_map>

using namespace std;


inline int getPrefix(string s0, string s1)
{
	int result = 0;
	for(int i = 0; i < s0.length() &&  i < s1.length(); i++)
	{
		if(s0[i] == s1[i]) result++;
		else break;
	}
	return result;
}

int main(){
	string s1 = "xiaomingxuhelloworld";
	string s2 = "xiaomingxuhelloeveryone";

	int prefix = getPrefix(s1, s2);

	cout << prefix << endl;
	string s3 = s1.substr(0, prefix);
	int index = s1.find('o');
	string s4 = s1.substr(index);
	cout << index << endl;
	cout << s3 << endl;
	cout << s4 << endl;
	cout << "================================================" << endl;

	unordered_map<string, int> myMap;
	myMap["a"] = 1;
	myMap["b"] = 2;
	myMap["c"] = 0;
	myMap.insert({"d", 0});
	myMap.insert({"a", 0});
	myMap["d"]++;

	for(auto x : myMap)
		cout << x.first << ": " << x.second << endl;
	cerr << "the size of myMap is: " << myMap.size() << endl;
	unordered_map<string, int>().swap(myMap);
	cerr << "the size of myMap is: " << myMap.size() << endl;

	return 0;
}

