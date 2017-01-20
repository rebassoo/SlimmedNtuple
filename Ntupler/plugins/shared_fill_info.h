struct FillInfo
{
	unsigned int fillNumber;
	unsigned int runMin, runMax;
	string alignmentTag;

	FillInfo(unsigned int _fn=0, unsigned int _rmin=0, unsigned int _rmax=0, const string &at="") :
		fillNumber(_fn), runMin(_rmin), runMax(_rmax), alignmentTag(at)
	{
	}
};

//----------------------------------------------------------------------------------------------------

struct FillInfoCollection : public vector<FillInfo> 
{
	FillInfo FindByRun(unsigned int run) const
	{
		for (const auto it : *this)
		{
			if (it.runMin <= run && it.runMax >= run)
				return it;
		}

		printf("ERROR in FillInfoCollection::FindByRun > run %u not found in the collection.\n", run);
		exit(1);

		return FillInfo();
	}
};

//----------------------------------------------------------------------------------------------------

FillInfoCollection fillInfoCollection;

//----------------------------------------------------------------------------------------------------

void InitFillInfoCollection()
{
	fillInfoCollection.push_back(FillInfo(4947, 273725, 273730, "margin/fill_4947/method x"));
	fillInfoCollection.push_back(FillInfo(4953, 274094, 274094, "margin/fill_4953/method x"));
	fillInfoCollection.push_back(FillInfo(4961, 274198, 274200, "margin/fill_4961/method x"));
	fillInfoCollection.push_back(FillInfo(4964, 274240, 274241, "margin/fill_4964/method x"));
	fillInfoCollection.push_back(FillInfo(4964, 274244, 274244, "no_margin/fill_4964/method x"));
	fillInfoCollection.push_back(FillInfo(4976, 274282, 274286, "margin/fill_4976/method x"));
	fillInfoCollection.push_back(FillInfo(4985, 274387, 274388, "no_margin/fill_4985/method x"));
	fillInfoCollection.push_back(FillInfo(4988, 274420, 274422, "no_margin/fill_4988/method x"));
	fillInfoCollection.push_back(FillInfo(4990, 274440, 274443, "no_margin/fill_4990/method x"));
	fillInfoCollection.push_back(FillInfo(5005, 274954, 274959, "no_margin/fill_5005/method x"));
	fillInfoCollection.push_back(FillInfo(5013, 274966, 274971, "no_margin/fill_5013/method x"));
	fillInfoCollection.push_back(FillInfo(5017, 274998, 275001, "no_margin/fill_5017/method x"));
	fillInfoCollection.push_back(FillInfo(5020, 275059, 275074, "no_margin/fill_5020/method x"));
	fillInfoCollection.push_back(FillInfo(5021, 275124, 275125, "no_margin/fill_5021/method x"));
	fillInfoCollection.push_back(FillInfo(5024, 275282, 275293, "no_margin/fill_5024/method x"));
	fillInfoCollection.push_back(FillInfo(5026, 275309, 275311, "no_margin/fill_5026/method x"));
	fillInfoCollection.push_back(FillInfo(5027, 275319, 275338, "no_margin/fill_5027/method x"));
	fillInfoCollection.push_back(FillInfo(5028, 275344, 275345, "no_margin/fill_5028/method x"));
	fillInfoCollection.push_back(FillInfo(5029, 275370, 275371, "no_margin/fill_5029/method x"));
	fillInfoCollection.push_back(FillInfo(5030, 275375, 275376, "no_margin/fill_5030/method x"));
	fillInfoCollection.push_back(FillInfo(5038, 275656, 275659, "no_margin/fill_5038/method x"));
	fillInfoCollection.push_back(FillInfo(5043, 275757, 275783, "no_margin/fill_5043/method x"));
	fillInfoCollection.push_back(FillInfo(5045, 275828, 275847, "no_margin/fill_5045/method x"));
}
