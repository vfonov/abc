// generated by Fast Light User Interface Designer (fluid) version 1.0109

#ifndef GUI_h
#define GUI_h
#include <FL/Fl.H>
#include "EMSParameters.h"
#include "EMSParametersXMLFile.h"
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Return_Button.H>
#include <FL/Fl_Tabs.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Output.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Browser.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Menu_Bar.H>

class GUI {
public:
  Fl_Double_Window* MakeWindow();
  Fl_Double_Window *aboutwindow;
private:
  void cb_Close_i(Fl_Button*, void*);
  static void cb_Close(Fl_Button*, void*);
public:
  Fl_Double_Window *orientwindow;
  Fl_Input *orientInput;
private:
  void cb_OK_i(Fl_Return_Button*, void*);
  static void cb_OK(Fl_Return_Button*, void*);
public:
  Fl_Double_Window *runwindow;
  Fl_Double_Window *mainwindow;
  Fl_Group *tab1;
private:
  void cb_Next_i(Fl_Button*, void*);
  static void cb_Next(Fl_Button*, void*);
public:
  Fl_Input *suffixInput;
  Fl_Output *atlasdirOutput;
  Fl_Input *prior1Input;
  Fl_Input *prior2Input;
  Fl_Input *prior3Input;
  Fl_Input *prior4Input;
private:
  void cb_Change_i(Fl_Button*, void*);
  static void cb_Change(Fl_Button*, void*);
public:
  Fl_Input *atlasOrientInput;
  Fl_Check_Button *warpAtlasButton;
  Fl_Input *warpInput_nx;
  Fl_Input *warpInput_ny;
  Fl_Input *warpInput_nz;
  Fl_Group *tab2;
  Fl_Browser *imageBrowser;
private:
  void cb_Add_i(Fl_Button*, void*);
  static void cb_Add(Fl_Button*, void*);
  void cb_Remove_i(Fl_Button*, void*);
  static void cb_Remove(Fl_Button*, void*);
public:
  Fl_Input *filterIterInput;
  Fl_Input *filterDtInput;
private:
  void cb_Next1_i(Fl_Button*, void*);
  static void cb_Next1(Fl_Button*, void*);
public:
  Fl_Input *biasDegreeInput;
  Fl_Browser *imageOrientBrowser;
private:
  void cb_Change1_i(Fl_Button*, void*);
  static void cb_Change1(Fl_Button*, void*);
public:
  Fl_Choice *filterMethodChoice;
  static Fl_Menu_Item menu_filterMethodChoice[];
  Fl_Group *tab3;
  Fl_Output *outdirOutput;
  Fl_Choice *formatChoice;
  static Fl_Menu_Item menu_formatChoice[];
private:
  void cb_Change2_i(Fl_Button*, void*);
  static void cb_Change2(Fl_Button*, void*);
  void cb_Next2_i(Fl_Button*, void*);
  static void cb_Next2(Fl_Button*, void*);
public:
  Fl_Group *tab4;
private:
  void cb_Run_i(Fl_Button*, void*);
  static void cb_Run(Fl_Button*, void*);
public:
  Fl_Check_Button *debugCheck;
  Fl_Check_Button *nowriteCheck;
  static Fl_Menu_Item menu_[];
private:
  void cb_New_i(Fl_Menu_*, void*);
  static void cb_New(Fl_Menu_*, void*);
  void cb_Load_i(Fl_Menu_*, void*);
  static void cb_Load(Fl_Menu_*, void*);
  void cb_Save_i(Fl_Menu_*, void*);
  static void cb_Save(Fl_Menu_*, void*);
  void cb_Quit_i(Fl_Menu_*, void*);
  static void cb_Quit(Fl_Menu_*, void*);
  void cb_About_i(Fl_Menu_*, void*);
  static void cb_About(Fl_Menu_*, void*);
public:
  void Initialize();
  EMSParameters::Pointer GetParameters();
  void SetGUIElements(EMSParameters::Pointer p);
};
#endif
