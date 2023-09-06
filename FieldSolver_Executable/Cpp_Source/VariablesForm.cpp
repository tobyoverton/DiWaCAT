#include "FieldFunctions.h"
#include "VariableForm.h"

using namespace System;
using namespace System::Windows::Forms;

[STAThread]

void UserInterface(array<String^>^ args) {
	Application::EnableVisualStyles();
	Application::SetCompatibleTextRenderingDefault(false);
	DiWaCAT_FieldSolverUI::VariablesForm form;
	Application::Run(%form);

	system("pause");

}
