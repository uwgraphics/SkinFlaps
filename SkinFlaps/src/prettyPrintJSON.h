#ifndef __PRETTY_PRINT_JSON__
#define __PRETTY_PRINT_JSON__

#include <string>

class prettyPrintJSON
{
public:
	prettyPrintJSON() {}
	~prettyPrintJSON() {}

	void convert(const char *packedJSONstring, std::string &prettyPrintOutput)
	{
		prettyPrintOutput.clear();
		std::string indent,*out = &prettyPrintOutput;
		indent.clear();
		const char *c = packedJSONstring;
		while (*c != '\0') {
			if (*c == '[') {
				out->push_back(*c);
				indent.push_back(' ');
				indent.push_back(' ');
				out->append("\n");
				out->append(indent);
			}
			else if (*c == ']') {
				indent.pop_back();
				indent.pop_back();
				out->append("\n");
				out->append(indent);
				out->push_back(*c);
			}
			else if (*c == '{') {
				out->push_back(*c);
				indent.push_back(' ');
				indent.push_back(' ');
				out->append("\n");
				out->append(indent);
			}
			else if (*c == '}') {
				indent.pop_back();
				indent.pop_back();
				out->append("\n");
				out->append(indent);
				out->push_back(*c);
			}
			else if (*c == ',') {
				out->push_back(*c);
				out->append("\n");
				out->append(indent);
			}
			else if (*c == ':') {
				out->push_back(*c);
//				out->push_back(' ');  // in Windows adds an unreadable space after colon in C:filename.  File still quite human readable.
			}
			else
				out->push_back(*c);
			++c;
		}
	};

};

#endif  // __PRETTY_PRINT_JSON__
