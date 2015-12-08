#include "IrrGui.h"
#include <stdlib.h>  // system, rand, srand, RAND_MAX

using namespace irr;
using namespace core;
using namespace scene;
using namespace video;
using namespace io;
using namespace gui;
using namespace chrono;

IrrGui::IrrGui(ChIrrApp* myapp, std::vector<Smarticle*> *mySmarticlesVec) {

	sv = mySmarticlesVec;
	// store pointer applicaiton
	app = myapp;
	// ..add a GUI slider to control friction


	if (bucketType == STRESSSTICK)
	{
		text_Angle = app->GetIGUIEnvironment()->addStaticText(L"Angle: 0, Increment: 0",
			rect<s32>(850, 45, 1050, 60), true);
	}
	else
	{
		text_Angle = app->GetIGUIEnvironment()->addStaticText(L"Angle: 0, Increment: 0",
			rect<s32>(850, 45, 1050, 60), true);
	}
		text_SmarticleAmt = app->GetIGUIEnvironment()->addStaticText(L"Layers: 0, Smarticles: 0",
			rect<s32>(850, 65, 1050, 80), true);

		text_Q = app->GetIGUIEnvironment()->addStaticText(L"Press Q to activate GUI1",
			rect<s32>(850, 85, 1050, 100), true);
		text_W = app->GetIGUIEnvironment()->addStaticText(L"Press W to activate GUI2",
			rect<s32>(850, 105, 1050, 120), true);
		text_E = app->GetIGUIEnvironment()->addStaticText(L"Press E to activate GUI3",
			rect<s32>(850, 125, 1050, 140), true);
		text_R = app->GetIGUIEnvironment()->addStaticText(L"Press R to vibrate around current angle",
			rect<s32>(850, 145, 1050, 160), true);
		text_T = app->GetIGUIEnvironment()->addStaticText(L"Press T to vibrate around specified angles",
			rect<s32>(850, 165, 1050, 180), true);

		angle1Input = app->GetIGUIEnvironment()->addEditBox(L"75",
			rect<s32>(1050, 165, 1100, 180), true);
		angle2Input = app->GetIGUIEnvironment()->addEditBox(L"75",
			rect<s32>(1100, 165, 1150, 180), true);

		text_Y = app->GetIGUIEnvironment()->addStaticText(L"Press Y to move cylinder away",
			rect<s32>(850, 185, 1050, 200), true);
		text_successful = app->GetIGUIEnvironment()->addStaticText(L"Successfully moved smarticles",
			rect<s32>(850, 205, 1050, 220), true);
		resetSuccessfulCount();
		
		

		inc = .1;
		rampAngle = 10 * CH_C_PI / 180;
		rampInc = 1.0 / 60.0;
		drum_freq = 1;
		drum_omega = drum_freq * 2 * CH_C_PI;
	}
	int IrrGui::successfulCount = 0;
	bool IrrGui::OnEvent(const SEvent& event) {
		// check if user moved the sliders with mouse..
		if (event.EventType == irr::EET_KEY_INPUT_EVENT && !event.KeyInput.PressedDown) {
			switch (event.KeyInput.Key)
			{
			case irr::KEY_KEY_Q:
				if (Smarticle::global_GUI_value != 1)
					Smarticle::global_GUI_value = 1;
				else
					Smarticle::global_GUI_value = 0;
				return true;
				break;

			case irr::KEY_KEY_W:
				if (Smarticle::global_GUI_value != 2)
					Smarticle::global_GUI_value = 2;
				else
					Smarticle::global_GUI_value = 0;
				return true;
				break;
			case irr::KEY_KEY_E:
				if (Smarticle::global_GUI_value != 3)
					Smarticle::global_GUI_value = 3;
				else
					Smarticle::global_GUI_value = 0;
				return true;
				break;
			case irr::KEY_KEY_R: //vibrate around current theta
				if (Smarticle::global_GUI_value != 4) //TODO create a boolean and then call a method later which performs this so we don't have to run through smarticle vec here!
				{

					double CurrTheta01;
					double CurrTheta12;
					Smarticle::global_GUI_value = 4;
					std::pair<double, double> angPair;
					for (size_t i = 0; i < sv->size(); i++) //get each particles current theta
					{

						Smarticle* sPtr = sv->at(i);

						MoveType currMoveType = sPtr->moveType;
						std::vector<std::pair<double, double>> *v;



						switch (currMoveType)
						{
						case 0:
							v = &sPtr->global;
							break;
						case 1:
							v = &sPtr->gui1;
							break;
						case 2:
							v = &sPtr->gui2;
							break;
						case 3:
							v = &sPtr->gui3;
							break;
						case 4:
							v = &sPtr->vib;
							break;
						case 5:
							v = &sPtr->ot;
							break;
						default:
							v = &sPtr->global;
							break;
						}


						CurrTheta01 = sPtr->GetAngle1();
						CurrTheta12 = sPtr->GetAngle2();
						//CurrTheta01 = v->at(sPtr->moveTypeIdxs.at(currMoveType)).first;
						//CurrTheta12 = v->at(sPtr->moveTypeIdxs.at(currMoveType)).second;
						sPtr->vib.clear();

						sPtr->vib.emplace_back(CurrTheta01, CurrTheta12);
						sPtr->vib.emplace_back(CurrTheta01 - vibAmp, CurrTheta12 - vibAmp);
						sPtr->vib.emplace_back(CurrTheta01, CurrTheta12);
						sPtr->vib.emplace_back(CurrTheta01 + vibAmp, CurrTheta12 + vibAmp);
					}
				}

				else
					Smarticle::global_GUI_value = 0;
				return true;
				break;

			case irr::KEY_KEY_T:
				if (Smarticle::global_GUI_value != 5)
				{
					std::pair<double, double> angPair;
					double ang1;
					double ang2;
					Smarticle::global_GUI_value = 5;
					for (size_t i = 0; i < sv->size(); i++) //get each particles current theta
					{
						Smarticle* sPtr = sv->at(i);

						ang1 = wcstod(angle1Input->getText(), NULL)*CH_C_PI / 180;
						ang2 = wcstod(angle2Input->getText(), NULL)*CH_C_PI / 180;
						sPtr->vib.clear();

						//in case strange values are written
						if (ang2 > CH_C_PI || ang2 < -CH_C_PI)
						{
							Smarticle::global_GUI_value = 0;
							return true;
							break;
						}
						if (ang2 > CH_C_PI || ang2 < -CH_C_PI)
						{
							Smarticle::global_GUI_value = 0;
							return true;
							break;
						}

						sPtr->vib.emplace_back(ang1, ang2);

						sPtr->vib.emplace_back(ang1 - vibAmp, ang2 - vibAmp);

						sPtr->vib.emplace_back(ang1, ang2);

						sPtr->vib.emplace_back(ang1 + vibAmp, ang2 + vibAmp);

					}
				}
				else
					Smarticle::global_GUI_value = 0;
				return true;
				break;

			case irr::KEY_KEY_Y: //remove container or floor
				if (bucket_exist)
				{
					switch (bucketType)
					{
					case CYLINDER: case STRESSSTICK: case HOOKRAISE: case KNOBCYLINDER:
						for (size_t i = 0; i < bucket_bod_vec.size(); i++)
						{
							bucket_bod_vec.at(i)->SetBodyFixed(false);

							bucket_bod_vec.at(i)->SetPos(ChVector<>(
								bucket_bod_vec.at(i)->GetPos().x,
								bucket_bod_vec.at(i)->GetPos().y + 1,
								bucket_bod_vec.at(i)->GetPos().z));
						}
						break;
					case HOPPER:
						bucket_bott->SetPos(ChVector<>(1, 0, 0));
						break;
					case RAMP:
						bucket_bott->SetPos(ChVector<>(1, 0, 0));
						break;
					}

					bucket_exist = false;
				}

				return true;
				break;
			case irr::KEY_KEY_J:
			{
				static bool ran = false;
				switch (bucketType)
				{
				case CYLINDER: case STRESSSTICK: case HOOKRAISE: case KNOBCYLINDER:
					for (size_t i = 0; i < bucket_bod_vec.size(); i++)
					{
						//bucket_bod_vec.at(i)->SetBodyFixed(false);
						bucket_bod_vec.at(i)->SetBodyFixed(true);

						double x = bucket_bod_vec.at(i)->GetPos().x;
						double y = bucket_bod_vec.at(i)->GetPos().y;
						double z = bucket_bod_vec.at(i)->GetPos().z;
						//GetLog() << "\n" << bucket_bod_vec.at(i).get_ptr()->GetPos() << "\n";

						double theta = atan2(y, x);
						//double theta = atan(x/y);
						double r = sqrt(x*x + y*y);
						bucket_bod_vec.at(i)->SetPos(ChVector<>(x - bucket_rad*.01*cos(theta), y - bucket_rad*.01*sin(theta), z));

					}

					break;
				case HOPPER:
					bucket_bott->SetPos(ChVector<>(100, 0, 0));
					break;
				case RAMP:
					bucket_bott->SetPos(ChVector<>(100, 0, 0));
					break;
				}

				//return true;
				break;

			}
			case irr::KEY_KEY_K:
			{
				static bool ran = false;
				switch (bucketType)
				{
				case CYLINDER: case STRESSSTICK: case HOOKRAISE: case KNOBCYLINDER:
					for (size_t i = 0; i < bucket_bod_vec.size(); i++)
					{
						//bucket_bod_vec.at(i)->SetBodyFixed(false);
						bucket_bod_vec.at(i)->SetBodyFixed(true);

						double x = bucket_bod_vec.at(i)->GetPos().x;
						double y = bucket_bod_vec.at(i)->GetPos().y;
						double z = bucket_bod_vec.at(i)->GetPos().z;
						//GetLog() << "\n" << bucket_bod_vec.at(i).get_ptr()->GetPos() << "\n";

						double theta = atan2(y, x);
						double r = sqrt(x*x + y*y);
						bucket_bod_vec.at(i)->SetPos(ChVector<>(x + bucket_rad*.01*cos(theta), y + bucket_rad*.01*sin(theta), z));

					}

					break;
				case HOPPER:
					bucket_bott->SetPos(ChVector<>(100, 0, 0));
					break;
				case RAMP:
					bucket_bott->SetPos(ChVector<>(100, 0, 0));
					break;
				}

				//return true;
				break;

			}
			case irr::KEY_KEY_1:
				switch (bucketType)
				{
				case DRUM:
					drum_freq = drum_freq - rampInc;
					break;
				case RAMP:

					rampAngle = Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()).x - rampInc * CH_C_PI / 180.0;
					bucket->SetRot(Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(
						Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()).x - rampInc * CH_C_PI / 180.0
						, 0, 0)));
					bucket_bott->SetRot(Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(
						Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()).x - rampInc * CH_C_PI / 180.0
						, 0, 0)));
					break;
				default:
					rampAngle = Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()).x - rampInc * CH_C_PI / 180.0;
					bucket->SetRot(Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(
						Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()).x - rampInc * CH_C_PI / 180.0
						, 0, 0)));
					break;
				}
				drawAngle();
				return true;
				break;
			case irr::KEY_KEY_2:			//increase angle of bucket by rampInc
				switch (bucketType)
				{
				case DRUM:
					drum_freq = drum_freq + rampInc;
					break;
				case RAMP:

					rampAngle = Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()).x + rampInc * CH_C_PI / 180.0;
					bucket->SetRot(Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(
						Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()).x + rampInc * CH_C_PI / 180.0
						, 0, 0)));
					bucket_bott->SetRot(Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(
						Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()).x + rampInc * CH_C_PI / 180.0
						, 0, 0)));
					break;
				default:
					rampAngle = Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()).x + rampInc * CH_C_PI / 180.0;
					bucket->SetRot(Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(
						Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()).x + rampInc * CH_C_PI / 180.0
						, 0, 0)));
					break;
				}
				drawAngle();
				return true;
				break;
			case irr::KEY_KEY_3:			//decrease rampInc
				rampInc = rampInc - 1.0 / 60.0;
				drawAngle();
				return true;
				break;
			case irr::KEY_KEY_4:			//increase rampInc
				rampInc = rampInc + 1.0 / 60.0;
				drawAngle();
				return true;
				break;
			case irr::KEY_KEY_C:			//decrease p
				for (size_t i = 0; i < sv->size(); i++) //get each particles current theta
				{
					Smarticle* sPtr = sv->at(i);
					sPtr->armsController->p_gain = sPtr->armsController->p_gain - inc;
				}
				drawAngle();
				return true;
				break;
			case irr::KEY_KEY_V:			//decrease i
				for (size_t i = 0; i < sv->size(); i++) //get each particles current theta
				{
					Smarticle* sPtr = sv->at(i);
					sPtr->armsController->i_gain = sPtr->armsController->i_gain - inc;
				}
				drawAngle();
				return true;
				break;
			case irr::KEY_KEY_B:			//decrease d
				for (size_t i = 0; i < sv->size(); i++) //get each particles current theta
				{
					Smarticle* sPtr = sv->at(i);
					sPtr->armsController->d_gain = sPtr->armsController->d_gain - inc;
				}
				drawAngle();
				return true;
				break;
			case irr::KEY_KEY_D:			//increase p
				for (size_t i = 0; i < sv->size(); i++) //get each particles current theta
				{
					Smarticle* sPtr = sv->at(i);
					sPtr->armsController->p_gain = sPtr->armsController->p_gain + inc;
				}
				drawAngle();
				return true;
				break;
			case irr::KEY_KEY_F:			//increase i
				for (size_t i = 0; i < sv->size(); i++) //get each particles current theta
				{
					Smarticle* sPtr = sv->at(i);
					sPtr->armsController->i_gain = sPtr->armsController->i_gain + inc;
				}
				drawAngle();
				return true;
				break;
			case irr::KEY_KEY_G:			//increase d
				for (size_t i = 0; i < sv->size(); i++) //get each particles current theta
				{
					Smarticle* sPtr = sv->at(i);
					sPtr->armsController->d_gain = sPtr->armsController->d_gain + inc;
				}
				drawAngle();
				return true;
				break;
			}

		}
		return false;
	}
	void IrrGui::drawSmarticleAmt(int numLayers)//nu
	{
		char message[100]; sprintf(message, "Layers: %d, Smarticles: %d, GUI: %d", numLayers, sv->size(), Smarticle::global_GUI_value);
		this->text_SmarticleAmt->setText(core::stringw(message).c_str());
	}
	void IrrGui::drawAngle()
	{
		if (bucketType == STRESSSTICK)
		{
			if (sv->size() > 0)
			{
				Smarticle* sPtr = sv->at(0);
				char message[100]; sprintf(message, "P:%g, I:%g, D:%g", sPtr->armsController->p_gain, sPtr->armsController->i_gain, sPtr->armsController->d_gain);
				this->text_Angle->setText(core::stringw(message).c_str());
			}
		}
		else if (bucketType == DRUM)
		{
			char message[100]; sprintf(message, "AngVel: %g rpm, Increment: %g", drum_freq * 60, rampInc * 60);
			this->text_Angle->setText(core::stringw(message).c_str());
		}
		else{
			char message[100]; sprintf(message, "Angle: %1.3g, Increment: %1.3g", Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()).x * 180 / CH_C_PI, rampInc);
			this->text_Angle->setText(core::stringw(message).c_str());
		}

	}
	void IrrGui::drawSuccessful()
	{
		int count = 0;
		for (size_t i = 0; i < sv->size(); i++) //get each particles current theta
		{
			if (sv->at(i)->successfulMotion) count++;
		}
		char message[100]; sprintf(message, "Successfully Moving: %d/%d", count, sv->size());
		this->text_successful->setText(core::stringw(message).c_str());
	}
	void IrrGui::drawSuccessful2()
	{
		char message[100]; sprintf(message, "Successfully Moving: %d/%d", successfulCount, sv->size());
		this->text_successful->setText(core::stringw(message).c_str());
	}
	void IrrGui::addSuccessful(Smarticle &sPtr)
	{
		if (sPtr.successfulMotion)
			successfulCount++;

	}

	void IrrGui::resetSuccessfulCount()
	{
		successfulCount = 0;
	}