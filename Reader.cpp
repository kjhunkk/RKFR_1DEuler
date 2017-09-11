#include "Reader.h"

Reader::Reader()
{
	_initial = _boundary = _timeInteg = "";
	_polyOrder = 0;
	_area = _sizeX = _CFL = _T = 0.0;
}

Reader::~Reader()
{

}

bool Reader::readFile(std::string name)
{
	std::ifstream file;
	file.open(name);
	bool file_open = file.is_open();

	// Clear stream state
	file.clear();

	// Read file
	while (file.is_open())
	{
		// Read text
		std::string text;
		std::getline(file, text);
		text.erase(std::remove(text.begin(), text.end(), ' '), text.end());

		// Read flux scheme
		if (text.find("$$FLUXSCHEME=", 0) != std::string::npos)
			_fluxScheme = text.substr(13);

		// Read correction function
		if (text.find("$$CORRECTIONFUNCTION=", 0) != std::string::npos)
			_correctFunc = text.substr(21);

		// Read limiter type
		if (text.find("$$LIMITER=", 0) != std::string::npos)
			_limiter = text.substr(10);

		// Read conserved variable number to limiting
		if (text.find("$$LIMITINGVARIABLE=", 0) != std::string::npos)
			_limitVar = std::stoi(text.substr(19));

		// Read initial condition
		if (text.find("$$INITIAL=", 0) != std::string::npos)
			_initial = text.substr(10);

		// Read boundary condition
		if (text.find("$$BOUNDARY=", 0) != std::string::npos)
			_boundary = text.substr(11);

		// Read time integration type
		if (text.find("$$TIMEINTEGRATION=", 0) != std::string::npos)
			_timeInteg = text.substr(18);

		// Read target time
		if (text.find("$$TARGETTIME=", 0) != std::string::npos)
			_T = std::stod(text.substr(13));

		// Read start point
		if (text.find("$$STARTPOINT=", 0) != std::string::npos)
			_startPoint = std::stod(text.substr(13));

		// Read area size
		if (text.find("$$AREA=", 0) != std::string::npos)
			_area = std::stod(text.substr(7));

		// Read grid size
		if (text.find("$$GRIDSIZE=", 0) != std::string::npos)
			_sizeX = std::stod(text.substr(11));

		// Read target time
		if (text.find("$$CFL=", 0) != std::string::npos)
			_CFL = std::stod(text.substr(6));

		// Read order of polynomial
		if (text.find("$$POLYNOMIALORDER=", 0) != std::string::npos)
			_polyOrder = std::stoi(text.substr(18));

		if (file.eof()) file.close();
	}

	// Print conditions
	print();

	return file_open;
}

void Reader::print() const
{
	std::cout << "----------Conditions----------\n";
	std::cout << "$$ Flux Scheme         : " << _fluxScheme << "\n";
	std::cout << "$$ Correction Function : " << _correctFunc << "\n";
	std::cout << "$$ Limiter             : " << _limiter << "\n";
	std::cout << "$$ Limiting variable   : " << _limitVar << "\n";
	std::cout << "$$ Initial condition   : " << _initial << "\n";
	std::cout << "$$ Boundary condition  : " << _boundary << "\n";
	std::cout << "$$ Time integration    : " << _timeInteg << "\n";
	std::cout << "$$ Start point         : " << _startPoint << "\n";
	std::cout << "$$ Area                : " << _area << "\n";
	std::cout << "$$ Grid size           : " << _sizeX << "\n";
	std::cout << "$$ Target time         : " << _T << "\n";
	std::cout << "$$ CFL number          : " << _CFL << "\n";
	std::cout << "$$ Order of polynomial : " << _polyOrder << "\n";
	std::cout << "------------------------------\n";
}