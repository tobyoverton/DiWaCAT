#include <msclr\marshal_cppstd.h>
#include <vector>
#include "FieldFunctions.h"
#pragma once

namespace DiWaCAT_FieldSolverUI {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;


	/// <summary>
	/// Summary for MyForm
	/// </summary>
	public ref class VariablesForm : public System::Windows::Forms::Form
	{
	public:
		VariablesForm(void)
		{
			InitializeComponent();
			//
			//TODO: Add the constructor code here
			//
		}

	protected:
		/// <summary>
		/// Clean up any resources being used.
		/// </summary>
		~VariablesForm()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::Button^  RunBox;
	protected:

	private: System::Windows::Forms::TextBox^  x0Box;
	private: System::Windows::Forms::TextBox^  y0Box;
	private: System::Windows::Forms::TextBox^  aBox;
	private: System::Windows::Forms::TextBox^  deltaBox;
	private: System::Windows::Forms::TextBox^  wBox;
	private: System::Windows::Forms::TextBox^  lBox;
	private: System::Windows::Forms::TextBox^  PermitivityBox;
	private: System::Windows::Forms::TextBox^  PermeabilityBox;








	private: System::Windows::Forms::Label^  label1;
	private: System::Windows::Forms::Label^  label2;
	private: System::Windows::Forms::Label^  label3;
	private: System::Windows::Forms::Label^  label4;
	private: System::Windows::Forms::Label^  label5;
	private: System::Windows::Forms::Label^  label6;
	private: System::Windows::Forms::Label^  label7;
	private: System::Windows::Forms::Label^  label8;
	private: System::Windows::Forms::Label^  label9;
	private: System::Windows::Forms::Label^  label10;
	private: System::Windows::Forms::TextBox^  nStepsBox;

	private: System::Windows::Forms::TextBox^  ModeAccuracyBox;

	private: System::Windows::Forms::TextBox^  RootBox;

	private: System::Windows::Forms::TextBox^  nYBox;

	private: System::Windows::Forms::TextBox^  nXBox;

	private: System::Windows::Forms::Label^  label11;
	private: System::Windows::Forms::Label^  label12;
	private: System::Windows::Forms::Label^  label13;
	private: System::Windows::Forms::Label^  label14;
	private: System::Windows::Forms::Label^  label15;
	private: System::Windows::Forms::Label^  label16;
	private: System::Windows::Forms::Label^  label17;
	private: System::Windows::Forms::Label^  label18;
	private: System::Windows::Forms::TextBox^  OutFileName;
	private: System::Windows::Forms::Label^  label19;
	private: System::Windows::Forms::Button^  DefaultParameters;
	private: System::Windows::Forms::Button^  CentreButton;
	private: System::Windows::Forms::Label^  label20;
	private: System::Windows::Forms::CheckBox^  ConvergenceTrueFalse;
	private: System::Windows::Forms::Label^  label21;
	private: System::Windows::Forms::ComboBox^  OrientationSelect;
	private: System::Windows::Forms::Button^  SelectFileButton;
	private: System::Windows::Forms::TextBox^  FileInName;
	private: System::Windows::Forms::TextBox^  OutFolderName;
	private: System::Windows::Forms::Button^  SelectOutFolderButton;
	private: System::Windows::Forms::OpenFileDialog^  InFileDialog;
	private: System::Windows::Forms::OpenFileDialog^  OutFolderDialog;
	private: System::Windows::Forms::Button^  QuitButton;
	private: System::Windows::Forms::Label^  label22;
	private: System::Windows::Forms::TextBox^  MaxZBox;
	private: System::Windows::Forms::ComboBox^  MaxZScale;
	private: System::Windows::Forms::Label^  CalcText;

	private: System::Windows::Forms::Button^  button1;
	private: System::Windows::Forms::Button^  DechirperButton;
	private: System::Windows::Forms::TabControl^  GeometryTab;
	private: System::Windows::Forms::TabPage^  PlanarTab;
	private: System::Windows::Forms::TabPage^  CircularTab;
	private: System::Windows::Forms::Label^  label23;
	private: System::Windows::Forms::TextBox^  x0BoxC;


	private: System::Windows::Forms::TextBox^  y0BoxC;
	private: System::Windows::Forms::TextBox^  aBoxC;
	private: System::Windows::Forms::TextBox^  deltaBoxC;
	private: System::Windows::Forms::ComboBox^  MaxZScaleC;
	private: System::Windows::Forms::TextBox^  MaxZBoxC;
	private: System::Windows::Forms::TextBox^  lBoxC;
	private: System::Windows::Forms::Label^  label24;
	private: System::Windows::Forms::TextBox^  PermitivityBoxC;
	private: System::Windows::Forms::TextBox^  PermeabilityBoxC;
	private: System::Windows::Forms::Label^  label25;
	private: System::Windows::Forms::Label^  label26;
	private: System::Windows::Forms::Label^  label27;
	private: System::Windows::Forms::Label^  label28;
	private: System::Windows::Forms::Label^  label31;
	private: System::Windows::Forms::CheckBox^  ConvergenceTrueFalseC;
	private: System::Windows::Forms::Label^  label32;
	private: System::Windows::Forms::Label^  label33;
	private: System::Windows::Forms::Label^  label34;
	private: System::Windows::Forms::Button^  CenterButtonC;
	private: System::Windows::Forms::Label^  label35;
	private: System::Windows::Forms::Button^  DefaultParametersC;
	private: System::Windows::Forms::TextBox^  nRBox;
	private: System::Windows::Forms::TextBox^  nTBox;
	private: System::Windows::Forms::TextBox^  RootBoxC;
	private: System::Windows::Forms::TextBox^  ModeAccuracyBoxC;
	private: System::Windows::Forms::TextBox^  nStepsBoxC;
	private: System::Windows::Forms::Label^  label36;
	private: System::Windows::Forms::Label^  label37;
	private: System::Windows::Forms::Label^  label38;
	private: System::Windows::Forms::Label^  label39;
	private: System::Windows::Forms::Label^  label40;
private: System::Windows::Forms::TabPage^  GreenTab;
private: System::Windows::Forms::CheckBox^  ConvergenceTrueFalseG;
private: System::Windows::Forms::Label^  label29;
private: System::Windows::Forms::TextBox^  x0BoxG;
private: System::Windows::Forms::Label^  label53;
private: System::Windows::Forms::TextBox^  y0BoxG;
private: System::Windows::Forms::Label^  label52;
private: System::Windows::Forms::TextBox^  aBoxG;
private: System::Windows::Forms::Label^  label51;
private: System::Windows::Forms::TextBox^  deltaBoxG;
private: System::Windows::Forms::Label^  label50;
private: System::Windows::Forms::ComboBox^  MaxZScaleG;
private: System::Windows::Forms::Label^  label49;
private: System::Windows::Forms::TextBox^  wBoxG;
private: System::Windows::Forms::TextBox^  nPointBoxG;
private: System::Windows::Forms::TextBox^  MaxZBoxG;
private: System::Windows::Forms::TextBox^  ModeAccBoxG;
private: System::Windows::Forms::Label^  label30;
private: System::Windows::Forms::TextBox^  RootBoxG;
private: System::Windows::Forms::TextBox^  PermitivityBoxG;
private: System::Windows::Forms::TextBox^  nYBoxG;
private: System::Windows::Forms::TextBox^  PermeabilityBoxG;
private: System::Windows::Forms::TextBox^  nXBoxG;
private: System::Windows::Forms::Label^  label41;
private: System::Windows::Forms::Label^  label48;
private: System::Windows::Forms::Label^  label42;
private: System::Windows::Forms::Label^  label47;
private: System::Windows::Forms::Label^  label43;
private: System::Windows::Forms::Label^  label46;
private: System::Windows::Forms::Label^  label44;
private: System::Windows::Forms::Label^  label45;











	protected:

	protected:

	private:
		/// <summary>
		/// Required designer variable.
		/// </summary>
		System::ComponentModel::Container ^components;

#pragma region Windows Form Designer generated code
		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		void InitializeComponent(void)
		{
			System::Windows::Forms::PictureBox^  pictureBox1;
			this->RunBox = (gcnew System::Windows::Forms::Button());
			this->x0Box = (gcnew System::Windows::Forms::TextBox());
			this->y0Box = (gcnew System::Windows::Forms::TextBox());
			this->aBox = (gcnew System::Windows::Forms::TextBox());
			this->deltaBox = (gcnew System::Windows::Forms::TextBox());
			this->wBox = (gcnew System::Windows::Forms::TextBox());
			this->lBox = (gcnew System::Windows::Forms::TextBox());
			this->PermitivityBox = (gcnew System::Windows::Forms::TextBox());
			this->PermeabilityBox = (gcnew System::Windows::Forms::TextBox());
			this->label1 = (gcnew System::Windows::Forms::Label());
			this->label2 = (gcnew System::Windows::Forms::Label());
			this->label3 = (gcnew System::Windows::Forms::Label());
			this->label4 = (gcnew System::Windows::Forms::Label());
			this->label5 = (gcnew System::Windows::Forms::Label());
			this->label6 = (gcnew System::Windows::Forms::Label());
			this->label7 = (gcnew System::Windows::Forms::Label());
			this->label8 = (gcnew System::Windows::Forms::Label());
			this->label9 = (gcnew System::Windows::Forms::Label());
			this->label10 = (gcnew System::Windows::Forms::Label());
			this->nStepsBox = (gcnew System::Windows::Forms::TextBox());
			this->ModeAccuracyBox = (gcnew System::Windows::Forms::TextBox());
			this->RootBox = (gcnew System::Windows::Forms::TextBox());
			this->nYBox = (gcnew System::Windows::Forms::TextBox());
			this->nXBox = (gcnew System::Windows::Forms::TextBox());
			this->label11 = (gcnew System::Windows::Forms::Label());
			this->label12 = (gcnew System::Windows::Forms::Label());
			this->label13 = (gcnew System::Windows::Forms::Label());
			this->label14 = (gcnew System::Windows::Forms::Label());
			this->label15 = (gcnew System::Windows::Forms::Label());
			this->label16 = (gcnew System::Windows::Forms::Label());
			this->label17 = (gcnew System::Windows::Forms::Label());
			this->label18 = (gcnew System::Windows::Forms::Label());
			this->OutFileName = (gcnew System::Windows::Forms::TextBox());
			this->label19 = (gcnew System::Windows::Forms::Label());
			this->DefaultParameters = (gcnew System::Windows::Forms::Button());
			this->CentreButton = (gcnew System::Windows::Forms::Button());
			this->label20 = (gcnew System::Windows::Forms::Label());
			this->ConvergenceTrueFalse = (gcnew System::Windows::Forms::CheckBox());
			this->label21 = (gcnew System::Windows::Forms::Label());
			this->OrientationSelect = (gcnew System::Windows::Forms::ComboBox());
			this->SelectFileButton = (gcnew System::Windows::Forms::Button());
			this->FileInName = (gcnew System::Windows::Forms::TextBox());
			this->OutFolderName = (gcnew System::Windows::Forms::TextBox());
			this->SelectOutFolderButton = (gcnew System::Windows::Forms::Button());
			this->InFileDialog = (gcnew System::Windows::Forms::OpenFileDialog());
			this->OutFolderDialog = (gcnew System::Windows::Forms::OpenFileDialog());
			this->QuitButton = (gcnew System::Windows::Forms::Button());
			this->label22 = (gcnew System::Windows::Forms::Label());
			this->MaxZBox = (gcnew System::Windows::Forms::TextBox());
			this->MaxZScale = (gcnew System::Windows::Forms::ComboBox());
			this->CalcText = (gcnew System::Windows::Forms::Label());
			this->button1 = (gcnew System::Windows::Forms::Button());
			this->DechirperButton = (gcnew System::Windows::Forms::Button());
			this->GeometryTab = (gcnew System::Windows::Forms::TabControl());
			this->PlanarTab = (gcnew System::Windows::Forms::TabPage());
			this->CircularTab = (gcnew System::Windows::Forms::TabPage());
			this->label23 = (gcnew System::Windows::Forms::Label());
			this->x0BoxC = (gcnew System::Windows::Forms::TextBox());
			this->y0BoxC = (gcnew System::Windows::Forms::TextBox());
			this->aBoxC = (gcnew System::Windows::Forms::TextBox());
			this->deltaBoxC = (gcnew System::Windows::Forms::TextBox());
			this->MaxZScaleC = (gcnew System::Windows::Forms::ComboBox());
			this->MaxZBoxC = (gcnew System::Windows::Forms::TextBox());
			this->lBoxC = (gcnew System::Windows::Forms::TextBox());
			this->label24 = (gcnew System::Windows::Forms::Label());
			this->PermitivityBoxC = (gcnew System::Windows::Forms::TextBox());
			this->PermeabilityBoxC = (gcnew System::Windows::Forms::TextBox());
			this->label25 = (gcnew System::Windows::Forms::Label());
			this->label26 = (gcnew System::Windows::Forms::Label());
			this->label27 = (gcnew System::Windows::Forms::Label());
			this->label28 = (gcnew System::Windows::Forms::Label());
			this->label31 = (gcnew System::Windows::Forms::Label());
			this->ConvergenceTrueFalseC = (gcnew System::Windows::Forms::CheckBox());
			this->label32 = (gcnew System::Windows::Forms::Label());
			this->label33 = (gcnew System::Windows::Forms::Label());
			this->label34 = (gcnew System::Windows::Forms::Label());
			this->CenterButtonC = (gcnew System::Windows::Forms::Button());
			this->label35 = (gcnew System::Windows::Forms::Label());
			this->DefaultParametersC = (gcnew System::Windows::Forms::Button());
			this->nRBox = (gcnew System::Windows::Forms::TextBox());
			this->nTBox = (gcnew System::Windows::Forms::TextBox());
			this->RootBoxC = (gcnew System::Windows::Forms::TextBox());
			this->ModeAccuracyBoxC = (gcnew System::Windows::Forms::TextBox());
			this->nStepsBoxC = (gcnew System::Windows::Forms::TextBox());
			this->label36 = (gcnew System::Windows::Forms::Label());
			this->label37 = (gcnew System::Windows::Forms::Label());
			this->label38 = (gcnew System::Windows::Forms::Label());
			this->label39 = (gcnew System::Windows::Forms::Label());
			this->label40 = (gcnew System::Windows::Forms::Label());
			this->GreenTab = (gcnew System::Windows::Forms::TabPage());
			this->label29 = (gcnew System::Windows::Forms::Label());
			this->x0BoxG = (gcnew System::Windows::Forms::TextBox());
			this->y0BoxG = (gcnew System::Windows::Forms::TextBox());
			this->aBoxG = (gcnew System::Windows::Forms::TextBox());
			this->deltaBoxG = (gcnew System::Windows::Forms::TextBox());
			this->MaxZScaleG = (gcnew System::Windows::Forms::ComboBox());
			this->wBoxG = (gcnew System::Windows::Forms::TextBox());
			this->MaxZBoxG = (gcnew System::Windows::Forms::TextBox());
			this->label30 = (gcnew System::Windows::Forms::Label());
			this->PermitivityBoxG = (gcnew System::Windows::Forms::TextBox());
			this->PermeabilityBoxG = (gcnew System::Windows::Forms::TextBox());
			this->label41 = (gcnew System::Windows::Forms::Label());
			this->label42 = (gcnew System::Windows::Forms::Label());
			this->label43 = (gcnew System::Windows::Forms::Label());
			this->label44 = (gcnew System::Windows::Forms::Label());
			this->label45 = (gcnew System::Windows::Forms::Label());
			this->label46 = (gcnew System::Windows::Forms::Label());
			this->label47 = (gcnew System::Windows::Forms::Label());
			this->label48 = (gcnew System::Windows::Forms::Label());
			this->nXBoxG = (gcnew System::Windows::Forms::TextBox());
			this->nYBoxG = (gcnew System::Windows::Forms::TextBox());
			this->RootBoxG = (gcnew System::Windows::Forms::TextBox());
			this->ModeAccBoxG = (gcnew System::Windows::Forms::TextBox());
			this->nPointBoxG = (gcnew System::Windows::Forms::TextBox());
			this->label49 = (gcnew System::Windows::Forms::Label());
			this->label50 = (gcnew System::Windows::Forms::Label());
			this->label51 = (gcnew System::Windows::Forms::Label());
			this->label52 = (gcnew System::Windows::Forms::Label());
			this->label53 = (gcnew System::Windows::Forms::Label());
			this->ConvergenceTrueFalseG = (gcnew System::Windows::Forms::CheckBox());
			pictureBox1 = (gcnew System::Windows::Forms::PictureBox());
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(pictureBox1))->BeginInit();
			this->GeometryTab->SuspendLayout();
			this->PlanarTab->SuspendLayout();
			this->CircularTab->SuspendLayout();
			this->GreenTab->SuspendLayout();
			this->SuspendLayout();
			// 
			// pictureBox1
			// 
			pictureBox1->Location = System::Drawing::Point(496, 153);
			pictureBox1->Name = L"pictureBox1";
			pictureBox1->Size = System::Drawing::Size(512, 218);
			pictureBox1->SizeMode = System::Windows::Forms::PictureBoxSizeMode::AutoSize;
			pictureBox1->TabIndex = 50;
			pictureBox1->TabStop = false;
			// 
			// RunBox
			// 
			this->RunBox->BackColor = System::Drawing::SystemColors::HotTrack;
			this->RunBox->Cursor = System::Windows::Forms::Cursors::Hand;
			this->RunBox->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 20));
			this->RunBox->ForeColor = System::Drawing::SystemColors::Control;
			this->RunBox->Location = System::Drawing::Point(821, 462);
			this->RunBox->Name = L"RunBox";
			this->RunBox->Size = System::Drawing::Size(75, 39);
			this->RunBox->TabIndex = 0;
			this->RunBox->Text = L"Run";
			this->RunBox->UseVisualStyleBackColor = false;
			this->RunBox->Click += gcnew System::EventHandler(this, &VariablesForm::RunButton_Click);
			// 
			// x0Box
			// 
			this->x0Box->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->x0Box->Location = System::Drawing::Point(141, 61);
			this->x0Box->Name = L"x0Box";
			this->x0Box->Size = System::Drawing::Size(87, 26);
			this->x0Box->TabIndex = 1;
			// 
			// y0Box
			// 
			this->y0Box->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->y0Box->Location = System::Drawing::Point(141, 93);
			this->y0Box->Name = L"y0Box";
			this->y0Box->Size = System::Drawing::Size(87, 26);
			this->y0Box->TabIndex = 2;
			// 
			// aBox
			// 
			this->aBox->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->aBox->Location = System::Drawing::Point(141, 125);
			this->aBox->Name = L"aBox";
			this->aBox->Size = System::Drawing::Size(87, 26);
			this->aBox->TabIndex = 3;
			// 
			// deltaBox
			// 
			this->deltaBox->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->deltaBox->Location = System::Drawing::Point(141, 157);
			this->deltaBox->Name = L"deltaBox";
			this->deltaBox->Size = System::Drawing::Size(87, 26);
			this->deltaBox->TabIndex = 4;
			// 
			// wBox
			// 
			this->wBox->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->wBox->Location = System::Drawing::Point(141, 189);
			this->wBox->Name = L"wBox";
			this->wBox->Size = System::Drawing::Size(87, 26);
			this->wBox->TabIndex = 5;
			// 
			// lBox
			// 
			this->lBox->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->lBox->Location = System::Drawing::Point(141, 221);
			this->lBox->Name = L"lBox";
			this->lBox->Size = System::Drawing::Size(87, 26);
			this->lBox->TabIndex = 6;
			// 
			// PermitivityBox
			// 
			this->PermitivityBox->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->PermitivityBox->Location = System::Drawing::Point(141, 253);
			this->PermitivityBox->Name = L"PermitivityBox";
			this->PermitivityBox->Size = System::Drawing::Size(87, 26);
			this->PermitivityBox->TabIndex = 7;
			// 
			// PermeabilityBox
			// 
			this->PermeabilityBox->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->PermeabilityBox->Location = System::Drawing::Point(141, 285);
			this->PermeabilityBox->Name = L"PermeabilityBox";
			this->PermeabilityBox->Size = System::Drawing::Size(87, 26);
			this->PermeabilityBox->TabIndex = 8;
			// 
			// label1
			// 
			this->label1->AutoSize = true;
			this->label1->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 17.5F));
			this->label1->Location = System::Drawing::Point(20, 25);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(195, 29);
			this->label1->TabIndex = 9;
			this->label1->Text = L"DLW Parameters";
			// 
			// label2
			// 
			this->label2->AutoSize = true;
			this->label2->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label2->Location = System::Drawing::Point(40, 64);
			this->label2->Name = L"label2";
			this->label2->Size = System::Drawing::Size(58, 20);
			this->label2->TabIndex = 10;
			this->label2->Text = L"x0 [cm]";
			// 
			// label3
			// 
			this->label3->AutoSize = true;
			this->label3->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label3->Location = System::Drawing::Point(40, 96);
			this->label3->Name = L"label3";
			this->label3->Size = System::Drawing::Size(58, 20);
			this->label3->TabIndex = 11;
			this->label3->Text = L"y0 [cm]";
			// 
			// label4
			// 
			this->label4->AutoSize = true;
			this->label4->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label4->Location = System::Drawing::Point(40, 128);
			this->label4->Name = L"label4";
			this->label4->Size = System::Drawing::Size(51, 20);
			this->label4->TabIndex = 12;
			this->label4->Text = L"a [cm]";
			// 
			// label5
			// 
			this->label5->AutoSize = true;
			this->label5->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label5->Location = System::Drawing::Point(40, 160);
			this->label5->Name = L"label5";
			this->label5->Size = System::Drawing::Size(51, 20);
			this->label5->TabIndex = 13;
			this->label5->Text = L"δ [cm]";
			// 
			// label6
			// 
			this->label6->AutoSize = true;
			this->label6->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label6->Location = System::Drawing::Point(40, 192);
			this->label6->Name = L"label6";
			this->label6->Size = System::Drawing::Size(53, 20);
			this->label6->TabIndex = 14;
			this->label6->Text = L"w [cm]";
			// 
			// label7
			// 
			this->label7->AutoSize = true;
			this->label7->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label7->Location = System::Drawing::Point(40, 224);
			this->label7->Name = L"label7";
			this->label7->Size = System::Drawing::Size(43, 20);
			this->label7->TabIndex = 15;
			this->label7->Text = L"L [m]";
			// 
			// label8
			// 
			this->label8->AutoSize = true;
			this->label8->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label8->Location = System::Drawing::Point(40, 256);
			this->label8->Name = L"label8";
			this->label8->Size = System::Drawing::Size(17, 20);
			this->label8->TabIndex = 16;
			this->label8->Text = L"ε";
			// 
			// label9
			// 
			this->label9->AutoSize = true;
			this->label9->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label9->Location = System::Drawing::Point(40, 288);
			this->label9->Name = L"label9";
			this->label9->Size = System::Drawing::Size(18, 20);
			this->label9->TabIndex = 17;
			this->label9->Text = L"μ";
			// 
			// label10
			// 
			this->label10->AutoSize = true;
			this->label10->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 17.5F));
			this->label10->Location = System::Drawing::Point(20, 361);
			this->label10->Name = L"label10";
			this->label10->Size = System::Drawing::Size(256, 29);
			this->label10->TabIndex = 18;
			this->label10->Text = L"Simulation Parameters";
			// 
			// nStepsBox
			// 
			this->nStepsBox->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->nStepsBox->Location = System::Drawing::Point(141, 554);
			this->nStepsBox->Name = L"nStepsBox";
			this->nStepsBox->Size = System::Drawing::Size(87, 26);
			this->nStepsBox->TabIndex = 23;
			// 
			// ModeAccuracyBox
			// 
			this->ModeAccuracyBox->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->ModeAccuracyBox->Location = System::Drawing::Point(141, 508);
			this->ModeAccuracyBox->Name = L"ModeAccuracyBox";
			this->ModeAccuracyBox->Size = System::Drawing::Size(87, 26);
			this->ModeAccuracyBox->TabIndex = 22;
			// 
			// RootBox
			// 
			this->RootBox->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->RootBox->Location = System::Drawing::Point(141, 466);
			this->RootBox->Name = L"RootBox";
			this->RootBox->Size = System::Drawing::Size(87, 26);
			this->RootBox->TabIndex = 21;
			// 
			// nYBox
			// 
			this->nYBox->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->nYBox->Location = System::Drawing::Point(141, 427);
			this->nYBox->Name = L"nYBox";
			this->nYBox->Size = System::Drawing::Size(87, 26);
			this->nYBox->TabIndex = 20;
			// 
			// nXBox
			// 
			this->nXBox->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->nXBox->Location = System::Drawing::Point(141, 393);
			this->nXBox->Name = L"nXBox";
			this->nXBox->Size = System::Drawing::Size(87, 26);
			this->nXBox->TabIndex = 19;
			// 
			// label11
			// 
			this->label11->AutoSize = true;
			this->label11->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label11->Location = System::Drawing::Point(40, 396);
			this->label11->Name = L"label11";
			this->label11->Size = System::Drawing::Size(29, 20);
			this->label11->TabIndex = 24;
			this->label11->Text = L"nX";
			// 
			// label12
			// 
			this->label12->AutoSize = true;
			this->label12->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label12->Location = System::Drawing::Point(40, 428);
			this->label12->Name = L"label12";
			this->label12->Size = System::Drawing::Size(29, 20);
			this->label12->TabIndex = 25;
			this->label12->Text = L"nY";
			// 
			// label13
			// 
			this->label13->AutoSize = true;
			this->label13->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label13->Location = System::Drawing::Point(40, 457);
			this->label13->Name = L"label13";
			this->label13->Size = System::Drawing::Size(73, 40);
			this->label13->TabIndex = 26;
			this->label13->Text = L"Root\r\nPrecision";
			// 
			// label14
			// 
			this->label14->AutoSize = true;
			this->label14->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label14->Location = System::Drawing::Point(40, 501);
			this->label14->Name = L"label14";
			this->label14->Size = System::Drawing::Size(100, 40);
			this->label14->TabIndex = 27;
			this->label14->Text = L"Mode\r\nAccuracy [%]";
			// 
			// label15
			// 
			this->label15->AutoSize = true;
			this->label15->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label15->Location = System::Drawing::Point(40, 557);
			this->label15->Name = L"label15";
			this->label15->Size = System::Drawing::Size(79, 20);
			this->label15->TabIndex = 28;
			this->label15->Text = L"No. Steps";
			// 
			// label16
			// 
			this->label16->AutoSize = true;
			this->label16->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label16->Location = System::Drawing::Point(16, 25);
			this->label16->Name = L"label16";
			this->label16->Size = System::Drawing::Size(102, 20);
			this->label16->TabIndex = 29;
			this->label16->Text = L"Beam In File:";
			// 
			// label17
			// 
			this->label17->AutoSize = true;
			this->label17->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label17->Location = System::Drawing::Point(443, 25);
			this->label17->Name = L"label17";
			this->label17->Size = System::Drawing::Size(134, 20);
			this->label17->TabIndex = 30;
			this->label17->Text = L"Beam Out Folder:";
			// 
			// label18
			// 
			this->label18->AutoSize = true;
			this->label18->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label18->Location = System::Drawing::Point(443, 51);
			this->label18->Name = L"label18";
			this->label18->Size = System::Drawing::Size(160, 20);
			this->label18->TabIndex = 31;
			this->label18->Text = L"Beam Out File Name:";
			// 
			// OutFileName
			// 
			this->OutFileName->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.25F));
			this->OutFileName->Location = System::Drawing::Point(609, 53);
			this->OutFileName->Name = L"OutFileName";
			this->OutFileName->Size = System::Drawing::Size(150, 20);
			this->OutFileName->TabIndex = 32;
			this->OutFileName->Text = L"DefaultBeamName";
			// 
			// label19
			// 
			this->label19->AutoSize = true;
			this->label19->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label19->Location = System::Drawing::Point(765, 53);
			this->label19->Name = L"label19";
			this->label19->Size = System::Drawing::Size(31, 20);
			this->label19->TabIndex = 33;
			this->label19->Text = L".h5";
			// 
			// DefaultParameters
			// 
			this->DefaultParameters->BackColor = System::Drawing::SystemColors::HotTrack;
			this->DefaultParameters->Cursor = System::Windows::Forms::Cursors::Hand;
			this->DefaultParameters->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->DefaultParameters->ForeColor = System::Drawing::SystemColors::Control;
			this->DefaultParameters->Location = System::Drawing::Point(253, 514);
			this->DefaultParameters->Name = L"DefaultParameters";
			this->DefaultParameters->Size = System::Drawing::Size(112, 52);
			this->DefaultParameters->TabIndex = 34;
			this->DefaultParameters->Text = L"Default\r\nParameters";
			this->DefaultParameters->UseVisualStyleBackColor = false;
			this->DefaultParameters->Click += gcnew System::EventHandler(this, &VariablesForm::DefaultParameters_Click);
			// 
			// CentreButton
			// 
			this->CentreButton->BackColor = System::Drawing::SystemColors::HotTrack;
			this->CentreButton->Cursor = System::Windows::Forms::Cursors::Hand;
			this->CentreButton->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->CentreButton->ForeColor = System::Drawing::SystemColors::Control;
			this->CentreButton->Location = System::Drawing::Point(253, 77);
			this->CentreButton->Name = L"CentreButton";
			this->CentreButton->Size = System::Drawing::Size(112, 29);
			this->CentreButton->TabIndex = 35;
			this->CentreButton->Text = L"Center Beam";
			this->CentreButton->UseVisualStyleBackColor = false;
			this->CentreButton->Click += gcnew System::EventHandler(this, &VariablesForm::CentreButton_Click);
			// 
			// label20
			// 
			this->label20->AutoSize = true;
			this->label20->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label20->Location = System::Drawing::Point(249, 393);
			this->label20->Name = L"label20";
			this->label20->Size = System::Drawing::Size(0, 20);
			this->label20->TabIndex = 36;
			// 
			// ConvergenceTrueFalse
			// 
			this->ConvergenceTrueFalse->AutoSize = true;
			this->ConvergenceTrueFalse->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->ConvergenceTrueFalse->Location = System::Drawing::Point(253, 410);
			this->ConvergenceTrueFalse->Name = L"ConvergenceTrueFalse";
			this->ConvergenceTrueFalse->Size = System::Drawing::Size(122, 24);
			this->ConvergenceTrueFalse->TabIndex = 38;
			this->ConvergenceTrueFalse->Text = L"Convergence";
			this->ConvergenceTrueFalse->UseVisualStyleBackColor = true;
			// 
			// label21
			// 
			this->label21->AutoSize = true;
			this->label21->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label21->Location = System::Drawing::Point(40, 320);
			this->label21->Name = L"label21";
			this->label21->Size = System::Drawing::Size(87, 20);
			this->label21->TabIndex = 39;
			this->label21->Text = L"Orientation";
			// 
			// OrientationSelect
			// 
			this->OrientationSelect->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 11));
			this->OrientationSelect->FormattingEnabled = true;
			this->OrientationSelect->Items->AddRange(gcnew cli::array< System::Object^  >(2) { L"Horizontal", L"Vertical" });
			this->OrientationSelect->Location = System::Drawing::Point(141, 317);
			this->OrientationSelect->Name = L"OrientationSelect";
			this->OrientationSelect->Size = System::Drawing::Size(174, 26);
			this->OrientationSelect->TabIndex = 40;
			this->OrientationSelect->Text = L"Horizontal";
			// 
			// SelectFileButton
			// 
			this->SelectFileButton->Location = System::Drawing::Point(376, 24);
			this->SelectFileButton->Name = L"SelectFileButton";
			this->SelectFileButton->Size = System::Drawing::Size(47, 25);
			this->SelectFileButton->TabIndex = 41;
			this->SelectFileButton->Text = L"Select";
			this->SelectFileButton->UseVisualStyleBackColor = true;
			this->SelectFileButton->Click += gcnew System::EventHandler(this, &VariablesForm::SelectFileButton_Click);
			// 
			// FileInName
			// 
			this->FileInName->Location = System::Drawing::Point(117, 25);
			this->FileInName->Name = L"FileInName";
			this->FileInName->Size = System::Drawing::Size(253, 20);
			this->FileInName->TabIndex = 42;
			this->FileInName->Text = L"No File Selected";
			// 
			// OutFolderName
			// 
			this->OutFolderName->Location = System::Drawing::Point(583, 27);
			this->OutFolderName->Name = L"OutFolderName";
			this->OutFolderName->Size = System::Drawing::Size(249, 20);
			this->OutFolderName->TabIndex = 43;
			this->OutFolderName->Text = L"No Folder Selected";
			// 
			// SelectOutFolderButton
			// 
			this->SelectOutFolderButton->Location = System::Drawing::Point(848, 24);
			this->SelectOutFolderButton->Name = L"SelectOutFolderButton";
			this->SelectOutFolderButton->Size = System::Drawing::Size(47, 25);
			this->SelectOutFolderButton->TabIndex = 44;
			this->SelectOutFolderButton->Text = L"Select";
			this->SelectOutFolderButton->UseVisualStyleBackColor = true;
			this->SelectOutFolderButton->Click += gcnew System::EventHandler(this, &VariablesForm::SelectOutFolderButton_Click);
			// 
			// InFileDialog
			// 
			this->InFileDialog->FileName = L"InFileDialog";
			// 
			// OutFolderDialog
			// 
			this->OutFolderDialog->AddExtension = false;
			this->OutFolderDialog->CheckFileExists = false;
			this->OutFolderDialog->FileName = L"Select Folder";
			// 
			// QuitButton
			// 
			this->QuitButton->BackColor = System::Drawing::SystemColors::ActiveCaptionText;
			this->QuitButton->Cursor = System::Windows::Forms::Cursors::Hand;
			this->QuitButton->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 20));
			this->QuitButton->ForeColor = System::Drawing::SystemColors::Control;
			this->QuitButton->Location = System::Drawing::Point(617, 462);
			this->QuitButton->Name = L"QuitButton";
			this->QuitButton->Size = System::Drawing::Size(75, 39);
			this->QuitButton->TabIndex = 45;
			this->QuitButton->Text = L"Quit";
			this->QuitButton->UseVisualStyleBackColor = false;
			this->QuitButton->Click += gcnew System::EventHandler(this, &VariablesForm::QuitButton_Click);
			// 
			// label22
			// 
			this->label22->AutoSize = true;
			this->label22->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label22->Location = System::Drawing::Point(40, 590);
			this->label22->Name = L"label22";
			this->label22->Size = System::Drawing::Size(128, 40);
			this->label22->TabIndex = 46;
			this->label22->Text = L"Max Longitudinal\r\nPosition";
			// 
			// MaxZBox
			// 
			this->MaxZBox->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->MaxZBox->Location = System::Drawing::Point(174, 590);
			this->MaxZBox->Name = L"MaxZBox";
			this->MaxZBox->Size = System::Drawing::Size(102, 26);
			this->MaxZBox->TabIndex = 47;
			// 
			// MaxZScale
			// 
			this->MaxZScale->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 11));
			this->MaxZScale->FormattingEnabled = true;
			this->MaxZScale->Items->AddRange(gcnew cli::array< System::Object^  >(2) { L"fs", L"m" });
			this->MaxZScale->Location = System::Drawing::Point(288, 590);
			this->MaxZScale->Name = L"MaxZScale";
			this->MaxZScale->Size = System::Drawing::Size(43, 26);
			this->MaxZScale->TabIndex = 48;
			this->MaxZScale->Text = L"m";
			// 
			// CalcText
			// 
			this->CalcText->AutoSize = true;
			this->CalcText->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->CalcText->Location = System::Drawing::Point(774, 509);
			this->CalcText->Name = L"CalcText";
			this->CalcText->Size = System::Drawing::Size(187, 20);
			this->CalcText->TabIndex = 49;
			this->CalcText->Text = L"Calculation Complete: No";
			// 
			// button1
			// 
			this->button1->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 30));
			this->button1->Location = System::Drawing::Point(658, 570);
			this->button1->Name = L"button1";
			this->button1->Size = System::Drawing::Size(221, 60);
			this->button1->TabIndex = 51;
			this->button1->Text = L"Run Test";
			this->button1->UseVisualStyleBackColor = true;
			this->button1->Click += gcnew System::EventHandler(this, &VariablesForm::button1_Click);
			// 
			// DechirperButton
			// 
			this->DechirperButton->BackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(128)),
				static_cast<System::Int32>(static_cast<System::Byte>(0)));
			this->DechirperButton->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->DechirperButton->Location = System::Drawing::Point(253, 174);
			this->DechirperButton->Name = L"DechirperButton";
			this->DechirperButton->Size = System::Drawing::Size(152, 28);
			this->DechirperButton->TabIndex = 52;
			this->DechirperButton->Text = L"CLARA Dechirper";
			this->DechirperButton->UseVisualStyleBackColor = false;
			this->DechirperButton->Click += gcnew System::EventHandler(this, &VariablesForm::DechirperButton_Click);
			// 
			// GeometryTab
			// 
			this->GeometryTab->Controls->Add(this->PlanarTab);
			this->GeometryTab->Controls->Add(this->CircularTab);
			this->GeometryTab->Controls->Add(this->GreenTab);
			this->GeometryTab->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->GeometryTab->Location = System::Drawing::Point(20, 84);
			this->GeometryTab->Name = L"GeometryTab";
			this->GeometryTab->SelectedIndex = 0;
			this->GeometryTab->Size = System::Drawing::Size(463, 676);
			this->GeometryTab->TabIndex = 53;
			// 
			// PlanarTab
			// 
			this->PlanarTab->Controls->Add(this->label1);
			this->PlanarTab->Controls->Add(this->DechirperButton);
			this->PlanarTab->Controls->Add(this->x0Box);
			this->PlanarTab->Controls->Add(this->y0Box);
			this->PlanarTab->Controls->Add(this->aBox);
			this->PlanarTab->Controls->Add(this->deltaBox);
			this->PlanarTab->Controls->Add(this->MaxZScale);
			this->PlanarTab->Controls->Add(this->wBox);
			this->PlanarTab->Controls->Add(this->MaxZBox);
			this->PlanarTab->Controls->Add(this->lBox);
			this->PlanarTab->Controls->Add(this->label22);
			this->PlanarTab->Controls->Add(this->PermitivityBox);
			this->PlanarTab->Controls->Add(this->PermeabilityBox);
			this->PlanarTab->Controls->Add(this->label2);
			this->PlanarTab->Controls->Add(this->label3);
			this->PlanarTab->Controls->Add(this->label4);
			this->PlanarTab->Controls->Add(this->label5);
			this->PlanarTab->Controls->Add(this->OrientationSelect);
			this->PlanarTab->Controls->Add(this->label6);
			this->PlanarTab->Controls->Add(this->label21);
			this->PlanarTab->Controls->Add(this->label7);
			this->PlanarTab->Controls->Add(this->ConvergenceTrueFalse);
			this->PlanarTab->Controls->Add(this->label8);
			this->PlanarTab->Controls->Add(this->label20);
			this->PlanarTab->Controls->Add(this->label9);
			this->PlanarTab->Controls->Add(this->CentreButton);
			this->PlanarTab->Controls->Add(this->label10);
			this->PlanarTab->Controls->Add(this->DefaultParameters);
			this->PlanarTab->Controls->Add(this->nXBox);
			this->PlanarTab->Controls->Add(this->nYBox);
			this->PlanarTab->Controls->Add(this->RootBox);
			this->PlanarTab->Controls->Add(this->ModeAccuracyBox);
			this->PlanarTab->Controls->Add(this->nStepsBox);
			this->PlanarTab->Controls->Add(this->label11);
			this->PlanarTab->Controls->Add(this->label15);
			this->PlanarTab->Controls->Add(this->label12);
			this->PlanarTab->Controls->Add(this->label14);
			this->PlanarTab->Controls->Add(this->label13);
			this->PlanarTab->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 15));
			this->PlanarTab->Location = System::Drawing::Point(4, 29);
			this->PlanarTab->Name = L"PlanarTab";
			this->PlanarTab->Padding = System::Windows::Forms::Padding(3);
			this->PlanarTab->Size = System::Drawing::Size(455, 643);
			this->PlanarTab->TabIndex = 0;
			this->PlanarTab->Text = L"Planar DLW";
			this->PlanarTab->UseVisualStyleBackColor = true;
			// 
			// CircularTab
			// 
			this->CircularTab->Controls->Add(this->label23);
			this->CircularTab->Controls->Add(this->x0BoxC);
			this->CircularTab->Controls->Add(this->y0BoxC);
			this->CircularTab->Controls->Add(this->aBoxC);
			this->CircularTab->Controls->Add(this->deltaBoxC);
			this->CircularTab->Controls->Add(this->MaxZScaleC);
			this->CircularTab->Controls->Add(this->MaxZBoxC);
			this->CircularTab->Controls->Add(this->lBoxC);
			this->CircularTab->Controls->Add(this->label24);
			this->CircularTab->Controls->Add(this->PermitivityBoxC);
			this->CircularTab->Controls->Add(this->PermeabilityBoxC);
			this->CircularTab->Controls->Add(this->label25);
			this->CircularTab->Controls->Add(this->label26);
			this->CircularTab->Controls->Add(this->label27);
			this->CircularTab->Controls->Add(this->label28);
			this->CircularTab->Controls->Add(this->label31);
			this->CircularTab->Controls->Add(this->ConvergenceTrueFalseC);
			this->CircularTab->Controls->Add(this->label32);
			this->CircularTab->Controls->Add(this->label33);
			this->CircularTab->Controls->Add(this->label34);
			this->CircularTab->Controls->Add(this->CenterButtonC);
			this->CircularTab->Controls->Add(this->label35);
			this->CircularTab->Controls->Add(this->DefaultParametersC);
			this->CircularTab->Controls->Add(this->nRBox);
			this->CircularTab->Controls->Add(this->nTBox);
			this->CircularTab->Controls->Add(this->RootBoxC);
			this->CircularTab->Controls->Add(this->ModeAccuracyBoxC);
			this->CircularTab->Controls->Add(this->nStepsBoxC);
			this->CircularTab->Controls->Add(this->label36);
			this->CircularTab->Controls->Add(this->label37);
			this->CircularTab->Controls->Add(this->label38);
			this->CircularTab->Controls->Add(this->label39);
			this->CircularTab->Controls->Add(this->label40);
			this->CircularTab->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 15));
			this->CircularTab->Location = System::Drawing::Point(4, 29);
			this->CircularTab->Name = L"CircularTab";
			this->CircularTab->Padding = System::Windows::Forms::Padding(3);
			this->CircularTab->Size = System::Drawing::Size(455, 643);
			this->CircularTab->TabIndex = 2;
			this->CircularTab->Text = L"Circular DLW";
			this->CircularTab->UseVisualStyleBackColor = true;
			// 
			// label23
			// 
			this->label23->AutoSize = true;
			this->label23->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 17.5F));
			this->label23->Location = System::Drawing::Point(20, 32);
			this->label23->Name = L"label23";
			this->label23->Size = System::Drawing::Size(195, 29);
			this->label23->TabIndex = 9;
			this->label23->Text = L"DLW Parameters";
			// 
			// x0BoxC
			// 
			this->x0BoxC->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->x0BoxC->Location = System::Drawing::Point(141, 70);
			this->x0BoxC->Name = L"x0BoxC";
			this->x0BoxC->Size = System::Drawing::Size(87, 26);
			this->x0BoxC->TabIndex = 1;
			// 
			// y0BoxC
			// 
			this->y0BoxC->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->y0BoxC->Location = System::Drawing::Point(141, 102);
			this->y0BoxC->Name = L"y0BoxC";
			this->y0BoxC->Size = System::Drawing::Size(87, 26);
			this->y0BoxC->TabIndex = 2;
			// 
			// aBoxC
			// 
			this->aBoxC->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->aBoxC->Location = System::Drawing::Point(141, 134);
			this->aBoxC->Name = L"aBoxC";
			this->aBoxC->Size = System::Drawing::Size(87, 26);
			this->aBoxC->TabIndex = 3;
			// 
			// deltaBoxC
			// 
			this->deltaBoxC->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->deltaBoxC->Location = System::Drawing::Point(141, 166);
			this->deltaBoxC->Name = L"deltaBoxC";
			this->deltaBoxC->Size = System::Drawing::Size(87, 26);
			this->deltaBoxC->TabIndex = 4;
			// 
			// MaxZScaleC
			// 
			this->MaxZScaleC->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 11));
			this->MaxZScaleC->FormattingEnabled = true;
			this->MaxZScaleC->Items->AddRange(gcnew cli::array< System::Object^  >(2) { L"fs", L"m" });
			this->MaxZScaleC->Location = System::Drawing::Point(288, 558);
			this->MaxZScaleC->Name = L"MaxZScaleC";
			this->MaxZScaleC->Size = System::Drawing::Size(43, 26);
			this->MaxZScaleC->TabIndex = 48;
			this->MaxZScaleC->Text = L"m";
			// 
			// MaxZBoxC
			// 
			this->MaxZBoxC->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->MaxZBoxC->Location = System::Drawing::Point(174, 558);
			this->MaxZBoxC->Name = L"MaxZBoxC";
			this->MaxZBoxC->Size = System::Drawing::Size(102, 26);
			this->MaxZBoxC->TabIndex = 47;
			// 
			// lBoxC
			// 
			this->lBoxC->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->lBoxC->Location = System::Drawing::Point(141, 203);
			this->lBoxC->Name = L"lBoxC";
			this->lBoxC->Size = System::Drawing::Size(87, 26);
			this->lBoxC->TabIndex = 6;
			// 
			// label24
			// 
			this->label24->AutoSize = true;
			this->label24->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label24->Location = System::Drawing::Point(40, 558);
			this->label24->Name = L"label24";
			this->label24->Size = System::Drawing::Size(128, 40);
			this->label24->TabIndex = 46;
			this->label24->Text = L"Max Longitudinal\r\nPosition";
			// 
			// PermitivityBoxC
			// 
			this->PermitivityBoxC->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->PermitivityBoxC->Location = System::Drawing::Point(141, 235);
			this->PermitivityBoxC->Name = L"PermitivityBoxC";
			this->PermitivityBoxC->Size = System::Drawing::Size(87, 26);
			this->PermitivityBoxC->TabIndex = 7;
			// 
			// PermeabilityBoxC
			// 
			this->PermeabilityBoxC->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->PermeabilityBoxC->Location = System::Drawing::Point(141, 267);
			this->PermeabilityBoxC->Name = L"PermeabilityBoxC";
			this->PermeabilityBoxC->Size = System::Drawing::Size(87, 26);
			this->PermeabilityBoxC->TabIndex = 8;
			// 
			// label25
			// 
			this->label25->AutoSize = true;
			this->label25->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label25->Location = System::Drawing::Point(40, 73);
			this->label25->Name = L"label25";
			this->label25->Size = System::Drawing::Size(58, 20);
			this->label25->TabIndex = 10;
			this->label25->Text = L"x0 [cm]";
			// 
			// label26
			// 
			this->label26->AutoSize = true;
			this->label26->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label26->Location = System::Drawing::Point(40, 105);
			this->label26->Name = L"label26";
			this->label26->Size = System::Drawing::Size(58, 20);
			this->label26->TabIndex = 11;
			this->label26->Text = L"y0 [cm]";
			// 
			// label27
			// 
			this->label27->AutoSize = true;
			this->label27->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label27->Location = System::Drawing::Point(40, 137);
			this->label27->Name = L"label27";
			this->label27->Size = System::Drawing::Size(51, 20);
			this->label27->TabIndex = 12;
			this->label27->Text = L"a [cm]";
			// 
			// label28
			// 
			this->label28->AutoSize = true;
			this->label28->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label28->Location = System::Drawing::Point(40, 169);
			this->label28->Name = L"label28";
			this->label28->Size = System::Drawing::Size(51, 20);
			this->label28->TabIndex = 13;
			this->label28->Text = L"δ [cm]";
			// 
			// label31
			// 
			this->label31->AutoSize = true;
			this->label31->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label31->Location = System::Drawing::Point(40, 206);
			this->label31->Name = L"label31";
			this->label31->Size = System::Drawing::Size(43, 20);
			this->label31->TabIndex = 15;
			this->label31->Text = L"L [m]";
			// 
			// ConvergenceTrueFalseC
			// 
			this->ConvergenceTrueFalseC->AutoSize = true;
			this->ConvergenceTrueFalseC->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->ConvergenceTrueFalseC->Location = System::Drawing::Point(253, 378);
			this->ConvergenceTrueFalseC->Name = L"ConvergenceTrueFalseC";
			this->ConvergenceTrueFalseC->Size = System::Drawing::Size(122, 24);
			this->ConvergenceTrueFalseC->TabIndex = 38;
			this->ConvergenceTrueFalseC->Text = L"Convergence";
			this->ConvergenceTrueFalseC->UseVisualStyleBackColor = true;
			// 
			// label32
			// 
			this->label32->AutoSize = true;
			this->label32->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label32->Location = System::Drawing::Point(40, 238);
			this->label32->Name = L"label32";
			this->label32->Size = System::Drawing::Size(17, 20);
			this->label32->TabIndex = 16;
			this->label32->Text = L"ε";
			// 
			// label33
			// 
			this->label33->AutoSize = true;
			this->label33->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label33->Location = System::Drawing::Point(249, 393);
			this->label33->Name = L"label33";
			this->label33->Size = System::Drawing::Size(0, 20);
			this->label33->TabIndex = 36;
			// 
			// label34
			// 
			this->label34->AutoSize = true;
			this->label34->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label34->Location = System::Drawing::Point(40, 270);
			this->label34->Name = L"label34";
			this->label34->Size = System::Drawing::Size(18, 20);
			this->label34->TabIndex = 17;
			this->label34->Text = L"μ";
			// 
			// CenterButtonC
			// 
			this->CenterButtonC->BackColor = System::Drawing::SystemColors::HotTrack;
			this->CenterButtonC->Cursor = System::Windows::Forms::Cursors::Hand;
			this->CenterButtonC->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->CenterButtonC->ForeColor = System::Drawing::SystemColors::Control;
			this->CenterButtonC->Location = System::Drawing::Point(253, 86);
			this->CenterButtonC->Name = L"CenterButtonC";
			this->CenterButtonC->Size = System::Drawing::Size(112, 29);
			this->CenterButtonC->TabIndex = 35;
			this->CenterButtonC->Text = L"Center Beam";
			this->CenterButtonC->UseVisualStyleBackColor = false;
			this->CenterButtonC->Click += gcnew System::EventHandler(this, &VariablesForm::CenterButtonC_Click);
			// 
			// label35
			// 
			this->label35->AutoSize = true;
			this->label35->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 17.5F));
			this->label35->Location = System::Drawing::Point(20, 325);
			this->label35->Name = L"label35";
			this->label35->Size = System::Drawing::Size(256, 29);
			this->label35->TabIndex = 18;
			this->label35->Text = L"Simulation Parameters";
			// 
			// DefaultParametersC
			// 
			this->DefaultParametersC->BackColor = System::Drawing::SystemColors::HotTrack;
			this->DefaultParametersC->Cursor = System::Windows::Forms::Cursors::Hand;
			this->DefaultParametersC->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->DefaultParametersC->ForeColor = System::Drawing::SystemColors::Control;
			this->DefaultParametersC->Location = System::Drawing::Point(253, 482);
			this->DefaultParametersC->Name = L"DefaultParametersC";
			this->DefaultParametersC->Size = System::Drawing::Size(112, 52);
			this->DefaultParametersC->TabIndex = 34;
			this->DefaultParametersC->Text = L"Default\r\nParameters";
			this->DefaultParametersC->UseVisualStyleBackColor = false;
			this->DefaultParametersC->Click += gcnew System::EventHandler(this, &VariablesForm::DefaultParametersC_Click);
			// 
			// nRBox
			// 
			this->nRBox->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->nRBox->Location = System::Drawing::Point(141, 361);
			this->nRBox->Name = L"nRBox";
			this->nRBox->Size = System::Drawing::Size(87, 26);
			this->nRBox->TabIndex = 19;
			// 
			// nTBox
			// 
			this->nTBox->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->nTBox->Location = System::Drawing::Point(141, 395);
			this->nTBox->Name = L"nTBox";
			this->nTBox->Size = System::Drawing::Size(87, 26);
			this->nTBox->TabIndex = 20;
			// 
			// RootBoxC
			// 
			this->RootBoxC->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->RootBoxC->Location = System::Drawing::Point(141, 434);
			this->RootBoxC->Name = L"RootBoxC";
			this->RootBoxC->Size = System::Drawing::Size(87, 26);
			this->RootBoxC->TabIndex = 21;
			// 
			// ModeAccuracyBoxC
			// 
			this->ModeAccuracyBoxC->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->ModeAccuracyBoxC->Location = System::Drawing::Point(141, 476);
			this->ModeAccuracyBoxC->Name = L"ModeAccuracyBoxC";
			this->ModeAccuracyBoxC->Size = System::Drawing::Size(87, 26);
			this->ModeAccuracyBoxC->TabIndex = 22;
			// 
			// nStepsBoxC
			// 
			this->nStepsBoxC->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->nStepsBoxC->Location = System::Drawing::Point(141, 522);
			this->nStepsBoxC->Name = L"nStepsBoxC";
			this->nStepsBoxC->Size = System::Drawing::Size(87, 26);
			this->nStepsBoxC->TabIndex = 23;
			// 
			// label36
			// 
			this->label36->AutoSize = true;
			this->label36->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label36->Location = System::Drawing::Point(40, 364);
			this->label36->Name = L"label36";
			this->label36->Size = System::Drawing::Size(63, 20);
			this->label36->TabIndex = 24;
			this->label36->Text = L"nRadial";
			// 
			// label37
			// 
			this->label37->AutoSize = true;
			this->label37->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label37->Location = System::Drawing::Point(40, 525);
			this->label37->Name = L"label37";
			this->label37->Size = System::Drawing::Size(79, 20);
			this->label37->TabIndex = 28;
			this->label37->Text = L"No. Steps";
			// 
			// label38
			// 
			this->label38->AutoSize = true;
			this->label38->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label38->Location = System::Drawing::Point(40, 396);
			this->label38->Name = L"label38";
			this->label38->Size = System::Drawing::Size(59, 20);
			this->label38->TabIndex = 25;
			this->label38->Text = L"nTheta";
			// 
			// label39
			// 
			this->label39->AutoSize = true;
			this->label39->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label39->Location = System::Drawing::Point(40, 469);
			this->label39->Name = L"label39";
			this->label39->Size = System::Drawing::Size(100, 40);
			this->label39->TabIndex = 27;
			this->label39->Text = L"Mode\r\nAccuracy [%]";
			// 
			// label40
			// 
			this->label40->AutoSize = true;
			this->label40->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label40->Location = System::Drawing::Point(40, 425);
			this->label40->Name = L"label40";
			this->label40->Size = System::Drawing::Size(73, 40);
			this->label40->TabIndex = 26;
			this->label40->Text = L"Root\r\nPrecision";
			// 
			// GreenTab
			// 
			this->GreenTab->Controls->Add(this->ConvergenceTrueFalseG);
			this->GreenTab->Controls->Add(this->label29);
			this->GreenTab->Controls->Add(this->x0BoxG);
			this->GreenTab->Controls->Add(this->label53);
			this->GreenTab->Controls->Add(this->y0BoxG);
			this->GreenTab->Controls->Add(this->label52);
			this->GreenTab->Controls->Add(this->aBoxG);
			this->GreenTab->Controls->Add(this->label51);
			this->GreenTab->Controls->Add(this->deltaBoxG);
			this->GreenTab->Controls->Add(this->label50);
			this->GreenTab->Controls->Add(this->MaxZScaleG);
			this->GreenTab->Controls->Add(this->label49);
			this->GreenTab->Controls->Add(this->wBoxG);
			this->GreenTab->Controls->Add(this->nPointBoxG);
			this->GreenTab->Controls->Add(this->MaxZBoxG);
			this->GreenTab->Controls->Add(this->ModeAccBoxG);
			this->GreenTab->Controls->Add(this->label30);
			this->GreenTab->Controls->Add(this->RootBoxG);
			this->GreenTab->Controls->Add(this->PermitivityBoxG);
			this->GreenTab->Controls->Add(this->nYBoxG);
			this->GreenTab->Controls->Add(this->PermeabilityBoxG);
			this->GreenTab->Controls->Add(this->nXBoxG);
			this->GreenTab->Controls->Add(this->label41);
			this->GreenTab->Controls->Add(this->label48);
			this->GreenTab->Controls->Add(this->label42);
			this->GreenTab->Controls->Add(this->label47);
			this->GreenTab->Controls->Add(this->label43);
			this->GreenTab->Controls->Add(this->label46);
			this->GreenTab->Controls->Add(this->label44);
			this->GreenTab->Controls->Add(this->label45);
			this->GreenTab->Location = System::Drawing::Point(4, 29);
			this->GreenTab->Name = L"GreenTab";
			this->GreenTab->Padding = System::Windows::Forms::Padding(3);
			this->GreenTab->Size = System::Drawing::Size(455, 643);
			this->GreenTab->TabIndex = 3;
			this->GreenTab->Text = L"1D Green\'s Function";
			this->GreenTab->UseVisualStyleBackColor = true;
			// 
			// label29
			// 
			this->label29->AutoSize = true;
			this->label29->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 17.5F));
			this->label29->Location = System::Drawing::Point(43, 49);
			this->label29->Name = L"label29";
			this->label29->Size = System::Drawing::Size(195, 29);
			this->label29->TabIndex = 61;
			this->label29->Text = L"DLW Parameters";
			// 
			// x0BoxG
			// 
			this->x0BoxG->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->x0BoxG->Location = System::Drawing::Point(164, 85);
			this->x0BoxG->Name = L"x0BoxG";
			this->x0BoxG->Size = System::Drawing::Size(87, 26);
			this->x0BoxG->TabIndex = 54;
			// 
			// y0BoxG
			// 
			this->y0BoxG->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->y0BoxG->Location = System::Drawing::Point(164, 117);
			this->y0BoxG->Name = L"y0BoxG";
			this->y0BoxG->Size = System::Drawing::Size(87, 26);
			this->y0BoxG->TabIndex = 55;
			// 
			// aBoxG
			// 
			this->aBoxG->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->aBoxG->Location = System::Drawing::Point(164, 149);
			this->aBoxG->Name = L"aBoxG";
			this->aBoxG->Size = System::Drawing::Size(87, 26);
			this->aBoxG->TabIndex = 56;
			// 
			// deltaBoxG
			// 
			this->deltaBoxG->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->deltaBoxG->Location = System::Drawing::Point(164, 181);
			this->deltaBoxG->Name = L"deltaBoxG";
			this->deltaBoxG->Size = System::Drawing::Size(87, 26);
			this->deltaBoxG->TabIndex = 57;
			// 
			// MaxZScaleG
			// 
			this->MaxZScaleG->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 11));
			this->MaxZScaleG->FormattingEnabled = true;
			this->MaxZScaleG->Items->AddRange(gcnew cli::array< System::Object^  >(2) { L"fs", L"m" });
			this->MaxZScaleG->Location = System::Drawing::Point(311, 546);
			this->MaxZScaleG->Name = L"MaxZScaleG";
			this->MaxZScaleG->Size = System::Drawing::Size(43, 26);
			this->MaxZScaleG->TabIndex = 82;
			this->MaxZScaleG->Text = L"m";
			// 
			// wBoxG
			// 
			this->wBoxG->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->wBoxG->Location = System::Drawing::Point(164, 213);
			this->wBoxG->Name = L"wBoxG";
			this->wBoxG->Size = System::Drawing::Size(87, 26);
			this->wBoxG->TabIndex = 58;
			// 
			// MaxZBoxG
			// 
			this->MaxZBoxG->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->MaxZBoxG->Location = System::Drawing::Point(197, 546);
			this->MaxZBoxG->Name = L"MaxZBoxG";
			this->MaxZBoxG->Size = System::Drawing::Size(102, 26);
			this->MaxZBoxG->TabIndex = 81;
			// 
			// label30
			// 
			this->label30->AutoSize = true;
			this->label30->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label30->Location = System::Drawing::Point(63, 546);
			this->label30->Name = L"label30";
			this->label30->Size = System::Drawing::Size(128, 40);
			this->label30->TabIndex = 80;
			this->label30->Text = L"Max Longitudinal\r\nPosition";
			// 
			// PermitivityBoxG
			// 
			this->PermitivityBoxG->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->PermitivityBoxG->Location = System::Drawing::Point(164, 245);
			this->PermitivityBoxG->Name = L"PermitivityBoxG";
			this->PermitivityBoxG->Size = System::Drawing::Size(87, 26);
			this->PermitivityBoxG->TabIndex = 59;
			// 
			// PermeabilityBoxG
			// 
			this->PermeabilityBoxG->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->PermeabilityBoxG->Location = System::Drawing::Point(164, 277);
			this->PermeabilityBoxG->Name = L"PermeabilityBoxG";
			this->PermeabilityBoxG->Size = System::Drawing::Size(87, 26);
			this->PermeabilityBoxG->TabIndex = 60;
			// 
			// label41
			// 
			this->label41->AutoSize = true;
			this->label41->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label41->Location = System::Drawing::Point(63, 88);
			this->label41->Name = L"label41";
			this->label41->Size = System::Drawing::Size(58, 20);
			this->label41->TabIndex = 62;
			this->label41->Text = L"x0 [cm]";
			// 
			// label42
			// 
			this->label42->AutoSize = true;
			this->label42->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label42->Location = System::Drawing::Point(63, 120);
			this->label42->Name = L"label42";
			this->label42->Size = System::Drawing::Size(58, 20);
			this->label42->TabIndex = 63;
			this->label42->Text = L"y0 [cm]";
			// 
			// label43
			// 
			this->label43->AutoSize = true;
			this->label43->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label43->Location = System::Drawing::Point(63, 152);
			this->label43->Name = L"label43";
			this->label43->Size = System::Drawing::Size(51, 20);
			this->label43->TabIndex = 64;
			this->label43->Text = L"a [cm]";
			// 
			// label44
			// 
			this->label44->AutoSize = true;
			this->label44->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label44->Location = System::Drawing::Point(63, 184);
			this->label44->Name = L"label44";
			this->label44->Size = System::Drawing::Size(51, 20);
			this->label44->TabIndex = 65;
			this->label44->Text = L"δ [cm]";
			// 
			// label45
			// 
			this->label45->AutoSize = true;
			this->label45->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label45->Location = System::Drawing::Point(63, 216);
			this->label45->Name = L"label45";
			this->label45->Size = System::Drawing::Size(53, 20);
			this->label45->TabIndex = 66;
			this->label45->Text = L"w [cm]";
			// 
			// label46
			// 
			this->label46->AutoSize = true;
			this->label46->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label46->Location = System::Drawing::Point(63, 248);
			this->label46->Name = L"label46";
			this->label46->Size = System::Drawing::Size(17, 20);
			this->label46->TabIndex = 67;
			this->label46->Text = L"ε";
			// 
			// label47
			// 
			this->label47->AutoSize = true;
			this->label47->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label47->Location = System::Drawing::Point(63, 280);
			this->label47->Name = L"label47";
			this->label47->Size = System::Drawing::Size(18, 20);
			this->label47->TabIndex = 68;
			this->label47->Text = L"μ";
			// 
			// label48
			// 
			this->label48->AutoSize = true;
			this->label48->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 17.5F));
			this->label48->Location = System::Drawing::Point(43, 317);
			this->label48->Name = L"label48";
			this->label48->Size = System::Drawing::Size(256, 29);
			this->label48->TabIndex = 69;
			this->label48->Text = L"Simulation Parameters";
			// 
			// nXBoxG
			// 
			this->nXBoxG->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->nXBoxG->Location = System::Drawing::Point(164, 349);
			this->nXBoxG->Name = L"nXBoxG";
			this->nXBoxG->Size = System::Drawing::Size(87, 26);
			this->nXBoxG->TabIndex = 70;
			// 
			// nYBoxG
			// 
			this->nYBoxG->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->nYBoxG->Location = System::Drawing::Point(164, 383);
			this->nYBoxG->Name = L"nYBoxG";
			this->nYBoxG->Size = System::Drawing::Size(87, 26);
			this->nYBoxG->TabIndex = 71;
			// 
			// RootBoxG
			// 
			this->RootBoxG->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->RootBoxG->Location = System::Drawing::Point(164, 422);
			this->RootBoxG->Name = L"RootBoxG";
			this->RootBoxG->Size = System::Drawing::Size(87, 26);
			this->RootBoxG->TabIndex = 72;
			// 
			// ModeAccBoxG
			// 
			this->ModeAccBoxG->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->ModeAccBoxG->Location = System::Drawing::Point(164, 464);
			this->ModeAccBoxG->Name = L"ModeAccBoxG";
			this->ModeAccBoxG->Size = System::Drawing::Size(87, 26);
			this->ModeAccBoxG->TabIndex = 73;
			// 
			// nPointBoxG
			// 
			this->nPointBoxG->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->nPointBoxG->Location = System::Drawing::Point(164, 510);
			this->nPointBoxG->Name = L"nPointBoxG";
			this->nPointBoxG->Size = System::Drawing::Size(87, 26);
			this->nPointBoxG->TabIndex = 74;
			// 
			// label49
			// 
			this->label49->AutoSize = true;
			this->label49->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label49->Location = System::Drawing::Point(63, 352);
			this->label49->Name = L"label49";
			this->label49->Size = System::Drawing::Size(29, 20);
			this->label49->TabIndex = 75;
			this->label49->Text = L"nX";
			// 
			// label50
			// 
			this->label50->AutoSize = true;
			this->label50->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label50->Location = System::Drawing::Point(63, 513);
			this->label50->Name = L"label50";
			this->label50->Size = System::Drawing::Size(81, 20);
			this->label50->TabIndex = 79;
			this->label50->Text = L"No. Points";
			// 
			// label51
			// 
			this->label51->AutoSize = true;
			this->label51->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label51->Location = System::Drawing::Point(63, 384);
			this->label51->Name = L"label51";
			this->label51->Size = System::Drawing::Size(29, 20);
			this->label51->TabIndex = 76;
			this->label51->Text = L"nY";
			// 
			// label52
			// 
			this->label52->AutoSize = true;
			this->label52->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label52->Location = System::Drawing::Point(63, 457);
			this->label52->Name = L"label52";
			this->label52->Size = System::Drawing::Size(100, 40);
			this->label52->TabIndex = 78;
			this->label52->Text = L"Mode\r\nAccuracy [%]";
			// 
			// label53
			// 
			this->label53->AutoSize = true;
			this->label53->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->label53->Location = System::Drawing::Point(63, 413);
			this->label53->Name = L"label53";
			this->label53->Size = System::Drawing::Size(73, 40);
			this->label53->TabIndex = 77;
			this->label53->Text = L"Root\r\nPrecision";
			// 
			// ConvergenceTrueFalseG
			// 
			this->ConvergenceTrueFalseG->AutoSize = true;
			this->ConvergenceTrueFalseG->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			this->ConvergenceTrueFalseG->Location = System::Drawing::Point(277, 364);
			this->ConvergenceTrueFalseG->Name = L"ConvergenceTrueFalseG";
			this->ConvergenceTrueFalseG->Size = System::Drawing::Size(122, 24);
			this->ConvergenceTrueFalseG->TabIndex = 54;
			this->ConvergenceTrueFalseG->Text = L"Convergence";
			this->ConvergenceTrueFalseG->UseVisualStyleBackColor = true;
			// 
			// VariablesForm
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(96, 96);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Dpi;
			this->ClientSize = System::Drawing::Size(1020, 772);
			this->Controls->Add(this->GeometryTab);
			this->Controls->Add(this->button1);
			this->Controls->Add(pictureBox1);
			this->Controls->Add(this->CalcText);
			this->Controls->Add(this->QuitButton);
			this->Controls->Add(this->SelectOutFolderButton);
			this->Controls->Add(this->OutFolderName);
			this->Controls->Add(this->FileInName);
			this->Controls->Add(this->SelectFileButton);
			this->Controls->Add(this->label19);
			this->Controls->Add(this->OutFileName);
			this->Controls->Add(this->label18);
			this->Controls->Add(this->label17);
			this->Controls->Add(this->label16);
			this->Controls->Add(this->RunBox);
			this->Name = L"VariablesForm";
			this->Text = L"DiWaCAT Simulation Parameters";
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(pictureBox1))->EndInit();
			this->GeometryTab->ResumeLayout(false);
			this->PlanarTab->ResumeLayout(false);
			this->PlanarTab->PerformLayout();
			this->CircularTab->ResumeLayout(false);
			this->CircularTab->PerformLayout();
			this->GreenTab->ResumeLayout(false);
			this->GreenTab->PerformLayout();
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion

	private: System::Void DefaultParameters_Click(System::Object^  sender, System::EventArgs^  e) {
		nStepsBox->Text = "0";
		nXBox->Text = "25";
		nYBox->Text = "25";
		RootBox->Text = "0.01";
		ModeAccuracyBox->Text = "0.5";
		MaxZBox->Text = "0";
		OrientationSelect->SelectedItem = "Horizontal";
		ConvergenceTrueFalse->Checked=true;
	}
	private: System::Void DefaultParametersC_Click(System::Object^  sender, System::EventArgs^  e) {
		nStepsBoxC->Text = "0";
		nRBox->Text = "10";
		nTBox->Text = "5";
		RootBoxC->Text = "0.1";
		ModeAccuracyBoxC->Text = "0.5";
		MaxZBoxC->Text = "0";
		ConvergenceTrueFalseC->Checked = true;
	}
	private: System::Void CentreButton_Click(System::Object^  sender, System::EventArgs^  e) {
		x0Box->Text = "0";
		y0Box->Text = "0";
	}
	private: System::Void CenterButtonC_Click(System::Object^  sender, System::EventArgs^  e) {
		x0BoxC->Text = "0";
		y0BoxC->Text = "0";
	}
	private: System::Void SelectFileButton_Click(System::Object^  sender, System::EventArgs^  e) {
		if (InFileDialog->ShowDialog() == System::Windows::Forms::DialogResult::OK) {
			FileInName->Text = InFileDialog->FileName;
		}
	}
	private: System::Void SelectOutFolderButton_Click(System::Object^  sender, System::EventArgs^  e) {
		if (OutFolderDialog->ShowDialog() == System::Windows::Forms::DialogResult::OK) {
			OutFolderName->Text = OutFolderDialog->FileName;
		}
	}

	private: System::Void RunButton_Click(System::Object^  sender, System::EventArgs^  e) {
		
		CalcText->Text = "Calculation Complete: No";
		std::vector<std::string> InputVariables;
		
		bool MissingValue = false;
		bool FileCorrect = true;

		//Check an input file is selected
		msclr::interop::marshal_context context;
		std::string InputFile = context.marshal_as<std::string>(FileInName->Text);
		if ((GeometryTab->SelectedTab != GeometryTab->TabPages["GreenTab"]) && ((InputFile.substr(InputFile.find_last_of(".") + 1) != "h5") || (InputFile == "No File Selected"))) {
			MessageBox::Show("Input File must be .h5", "Input File Error", MessageBoxButtons::OK, MessageBoxIcon::Error);
			FileCorrect = false;
		}

		//Set up the output file - make sure the folder is a folder and not a file and make sure the file has no extension
		std::string OutputFolder = context.marshal_as<std::string>(OutFolderName->Text);
		std::string OutputFileName = context.marshal_as<std::string>(OutFileName->Text);
		OutputFolder = OutputFolder.substr(0, OutputFolder.find_last_of('\\'));
		OutputFileName = OutputFileName.substr(0, OutputFileName.find_first_of('.'));
		std::string OutputFile = OutputFolder + "\\" + OutputFileName + ".h5";

		//Check if any boxes are empty (messy coding but it works)
		if (GeometryTab->SelectedTab == GeometryTab->TabPages["PlanarTab"]) {
			if ((x0Box->Text == "") || (y0Box->Text == "") || (PermitivityBox->Text == "") || (PermeabilityBox->Text == "") || (aBox->Text == "") || (deltaBox->Text == "") || (wBox->Text == "") || (nXBox->Text == "") || (nYBox->Text == "") || (RootBox->Text == "") || (ModeAccuracyBox->Text == "") || (lBox->Text == "") || (nStepsBox->Text == "") || (MaxZBox->Text == "")) {
				MessageBox::Show("Missing Value(s)", "Empty Value", MessageBoxButtons::OK, MessageBoxIcon::Error);
				MissingValue = true;
			}
		}
		if (GeometryTab->SelectedTab == GeometryTab->TabPages["CircularTab"]) {
			if ((x0BoxC->Text == "") || (y0BoxC->Text == "") || (PermitivityBoxC->Text == "") || (PermeabilityBoxC->Text == "") || (aBoxC->Text == "") || (deltaBoxC->Text == "") || (nRBox->Text == "") || (nTBox->Text == "") || (RootBoxC->Text == "") || (ModeAccuracyBoxC->Text == "") || (lBoxC->Text == "") || (nStepsBoxC->Text == "") || (MaxZBoxC->Text == "")) {
				MessageBox::Show("Missing Value(s)", "Empty Value", MessageBoxButtons::OK, MessageBoxIcon::Error);
				MissingValue = true;
			}
		}
		if (GeometryTab->SelectedTab == GeometryTab->TabPages["GreenTab"]) {
			if ((x0BoxG->Text == "") || (y0BoxG->Text == "") || (PermitivityBoxG->Text == "") || (PermeabilityBoxG->Text == "") || (aBoxG->Text == "") || (deltaBoxG->Text == "") || (wBoxG->Text=="") || (nXBoxG->Text == "") || (nYBoxG->Text == "") || (RootBoxG->Text == "") ||(ModeAccBoxG->Text == "")|| nPointBoxG->Text == "" || MaxZBoxG->Text == "") {
				MessageBox::Show("Missing Value(s)", "Empty Value", MessageBoxButtons::OK, MessageBoxIcon::Error);
				MissingValue = true;
			}
		}

		if(MissingValue==false && FileCorrect==true) {
			//Add the variables in the order needed for DiWaCAT
			if (GeometryTab->SelectedTab == GeometryTab->TabPages["PlanarTab"]) {
				if (OrientationSelect->Text == "Vertical") {
					InputVariables.push_back("v");
				}
				else {
					InputVariables.push_back("h");
				}
				InputVariables.push_back(context.marshal_as<std::string>(x0Box->Text));
				InputVariables.push_back(context.marshal_as<std::string>(y0Box->Text));
				InputVariables.push_back(context.marshal_as<std::string>(PermitivityBox->Text));
				InputVariables.push_back(context.marshal_as<std::string>(PermeabilityBox->Text));
				InputVariables.push_back(context.marshal_as<std::string>(aBox->Text));
				InputVariables.push_back(context.marshal_as<std::string>(deltaBox->Text));
				InputVariables.push_back(context.marshal_as<std::string>(wBox->Text));
				InputVariables.push_back(context.marshal_as<std::string>(nXBox->Text));
				InputVariables.push_back(context.marshal_as<std::string>(nYBox->Text));
				InputVariables.push_back(context.marshal_as<std::string>(RootBox->Text));
				InputVariables.push_back(context.marshal_as<std::string>(ModeAccuracyBox->Text));
				InputVariables.push_back(context.marshal_as<std::string>(lBox->Text));
				InputVariables.push_back(context.marshal_as<std::string>(nStepsBox->Text));
				//If longitudinal position set to seconds, convert to cm for DiWaCAT
				if (MaxZScale->Text == "fs") {
					double MaxZ = std::stold(context.marshal_as<std::string>(MaxZBox->Text))*2.99792e8*1e-15;
					InputVariables.push_back(std::to_string(MaxZ));
				}
				else {
					InputVariables.push_back(context.marshal_as<std::string>(MaxZBox->Text));
				}
				if (ConvergenceTrueFalse->Checked) {
					InputVariables.push_back("t");
				}
				else {
					InputVariables.push_back("f");
				}
				//Calculate field for the variables given in UI and file set
				std::cout << "Calling Field Solver Function" << std::endl;
				DiWaCATField_Planar(InputVariables, InputFile, OutputFile);
				CalcText->Text = "Calculation Complete: Yes";
			}
			if (GeometryTab->SelectedTab == GeometryTab->TabPages["CircularTab"]) {
				InputVariables.push_back(context.marshal_as<std::string>(x0BoxC->Text));
				InputVariables.push_back(context.marshal_as<std::string>(y0BoxC->Text));
				InputVariables.push_back(context.marshal_as<std::string>(PermitivityBoxC->Text));
				InputVariables.push_back(context.marshal_as<std::string>(PermeabilityBoxC->Text));
				InputVariables.push_back(context.marshal_as<std::string>(aBoxC->Text));
				InputVariables.push_back(context.marshal_as<std::string>(deltaBoxC->Text));
				InputVariables.push_back(context.marshal_as<std::string>(nRBox->Text));
				InputVariables.push_back(context.marshal_as<std::string>(nTBox->Text));
				InputVariables.push_back(context.marshal_as<std::string>(RootBoxC->Text));
				InputVariables.push_back(context.marshal_as<std::string>(ModeAccuracyBoxC->Text));
				InputVariables.push_back(context.marshal_as<std::string>(lBoxC->Text));
				InputVariables.push_back(context.marshal_as<std::string>(nStepsBoxC->Text));
				//If longitudinal position set to seconds, convert to cm for DiWaCAT
				if (MaxZScaleC->Text == "fs") {
					double MaxZ = std::stold(context.marshal_as<std::string>(MaxZBoxC->Text))*2.99792e8*1e-15;
					InputVariables.push_back(std::to_string(MaxZ));
				}
				else {
					InputVariables.push_back(context.marshal_as<std::string>(MaxZBoxC->Text));
				}
				if (ConvergenceTrueFalseC->Checked) {
					InputVariables.push_back("t");
				}
				else {
					InputVariables.push_back("f");
				}
				//Calculate field for the variables given in UI and file set
				std::cout << "Calling Field Solver Function" << std::endl;
				DiWaCATField_Circ(InputVariables, InputFile, OutputFile);
				CalcText->Text = "Calculation Complete: Yes";
			}
			if (GeometryTab->SelectedTab == GeometryTab->TabPages["GreenTab"]) {
				InputVariables.push_back(context.marshal_as<std::string>(x0BoxG->Text));
				InputVariables.push_back(context.marshal_as<std::string>(y0BoxG->Text));
				InputVariables.push_back(context.marshal_as<std::string>(PermitivityBoxG->Text));
				InputVariables.push_back(context.marshal_as<std::string>(PermeabilityBoxG->Text));
				InputVariables.push_back(context.marshal_as<std::string>(aBoxG->Text));
				InputVariables.push_back(context.marshal_as<std::string>(deltaBoxG->Text));
				InputVariables.push_back(context.marshal_as<std::string>(wBoxG->Text));
				InputVariables.push_back(context.marshal_as<std::string>(nXBoxG->Text));
				InputVariables.push_back(context.marshal_as<std::string>(nYBoxG->Text));
				InputVariables.push_back(context.marshal_as<std::string>(RootBoxG->Text));
				InputVariables.push_back(context.marshal_as<std::string>(ModeAccBoxG->Text));
				InputVariables.push_back(context.marshal_as<std::string>(nPointBoxG->Text));
				if (MaxZScaleG->Text == "fs") {
					double MaxZ = std::stold(context.marshal_as<std::string>(MaxZBoxG->Text))*2.99792e10*1e-15;
					InputVariables.push_back(std::to_string(MaxZ));
				}
				else {
					double MaxZ = std::stold(context.marshal_as<std::string>(MaxZBoxG->Text))*1e2;
					InputVariables.push_back(std::to_string(MaxZ));
				}
				if (ConvergenceTrueFalseG->Checked) {
					InputVariables.push_back("t");
				}
				else {
					InputVariables.push_back("f");
				}
				//Calculate field for the variables given in UI and file set
				std::cout << "Calling 1D Green's Function Calculator" << std::endl;
				DiWaCATField_1DGreen_Planar(InputVariables, OutputFile);
				CalcText->Text = "Calculation Complete: Yes";
			}
		}
	}
	private: System::Void QuitButton_Click(System::Object^  sender, System::EventArgs^  e) {
		Application::Exit();
	}

	private: System::Void DechirperButton_Click(System::Object^  sender, System::EventArgs^  e) {
		PermitivityBox->Text = "4";
		PermeabilityBox->Text = "1";
		deltaBox->Text = "0.02";
		wBox->Text = "2";
		lBox->Text = "0.2";
	}

	private: System::Void button1_Click(System::Object^  sender, System::EventArgs^  e) {
		CentreButton->PerformClick();
		PermitivityBox->Text = "4";
		PermeabilityBox->Text = "1";
		aBox->Text = "0.1";
		deltaBox->Text = "0.1";
		wBox->Text = "1";
		nXBox->Text = "2";
		nYBox->Text = "2";
		RootBox->Text = "0.01";
		ModeAccuracyBox->Text = "0.01";
		lBox->Text = "1";
		nStepsBox->Text = "0";
		MaxZBox->Text = "0";
		OutFolderName->Text = ".\\";
		OutFileName->Text = "TestBeamOut";
		FileInName->Text = ".\\FieldSolver_Executable\\600fs_Gaus.h5";
		RunBox->PerformClick();
	}

};
}
