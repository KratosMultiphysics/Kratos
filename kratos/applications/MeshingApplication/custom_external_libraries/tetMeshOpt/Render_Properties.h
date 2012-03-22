#pragma once

namespace PrSpain {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;

	/// <summary>
	/// Summary for Render_Properties
	/// </summary>
	public ref class Render_Properties : public System::Windows::Forms::Form
	{
	public:
		Render_Properties(void)
		{
			InitializeComponent();
			//
			//TODO: Add the constructor code here
			//
		}
		void setMesh()
		{
			//---------------
		}
	protected:
		/// <summary>
		/// Clean up any resources being used.
		/// </summary>
		~Render_Properties()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::Label^  label1;
	protected: 
	private: System::Windows::Forms::NumericUpDown^  numericUpDown1;

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
			this->label1 = (gcnew System::Windows::Forms::Label());
			this->numericUpDown1 = (gcnew System::Windows::Forms::NumericUpDown());
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->numericUpDown1))->BeginInit();
			this->SuspendLayout();
			// 
			// label1
			// 
			this->label1->AutoSize = true;
			this->label1->Location = System::Drawing::Point(13, 14);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(75, 13);
			this->label1->TabIndex = 3;
			this->label1->Text = L"Element Scale";
			// 
			// numericUpDown1
			// 
			this->numericUpDown1->Location = System::Drawing::Point(96, 12);
			this->numericUpDown1->Name = L"numericUpDown1";
			this->numericUpDown1->Size = System::Drawing::Size(78, 20);
			this->numericUpDown1->TabIndex = 2;
			this->numericUpDown1->Click += gcnew System::EventHandler(this, &Render_Properties::numericUpDown1_Click);
			// 
			// Render_Properties
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(197, 262);
			this->Controls->Add(this->label1);
			this->Controls->Add(this->numericUpDown1);
			this->Name = L"Render_Properties";
			this->Text = L"Render_Properties";
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->numericUpDown1))->EndInit();
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion
	private: System::Void numericUpDown1_Click(System::Object^  sender, System::EventArgs^  e) 
			 {
				 //----

			 }
	};
}
