using System;
using System.Collections.Generic;

//--------------------- Copyright Block ----------------------
/* 

PrayTime.cs: Prayer Times Calculator
Copyright (C) 2013
By: Ali Alghamdi
Original JS Code By: Hamid Zarrabi-Zadeh
praytimes.org

*/

namespace newpackage
{



	public class PrayTime
	{

		// ---------------------- Global Variables --------------------
		private int calcMethod; // caculation method
		private int asrJuristic; // Juristic method for Asr
		private int dhuhrMinutes; // minutes after mid-day for Dhuhr
		private int adjustHighLats; // adjusting method for higher latitudes
		private int timeFormat; // time format
		private double lat; // latitude
		private double lng; // longitude
		private double timeZone; // time-zone
		private double JDate_Renamed; // Julian date
		// ------------------------------------------------------------
		// Calculation Methods
		private int Jafari_Renamed; // Ithna Ashari
		private int Karachi_Renamed; // University of Islamic Sciences, Karachi
		private int ISNA_Renamed; // Islamic Society of North America (ISNA)
		private int MWL_Renamed; // Muslim World League (MWL)
		private int Makkah_Renamed; // Umm al-Qura, Makkah
		private int Egypt_Renamed; // Egyptian General Authority of Survey
		private int Custom_Renamed; // Custom Setting
		private int Tehran_Renamed; // Institute of Geophysics, University of Tehran
		// Juristic Methods
		private int Shafii_Renamed; // Shafii (standard)
		private int Hanafi_Renamed; // Hanafi
		// Adjusting Methods for Higher Latitudes
		private int None_Renamed; // No adjustment
		private int MidNight_Renamed; // middle of night
		private int OneSeventh_Renamed; // 1/7th of night
		private int AngleBased_Renamed; // angle/60th of night
		// Time Formats
		private int Time24_Renamed; // 24-hour format
		private int Time12_Renamed; // 12-hour format
		private int Time12NS_Renamed; // 12-hour format with no suffix
		private int Floating_Renamed; // floating point number
		// Time Names
		private List<string> timeNames;
		private string InvalidTime; // The string used for invalid times
		// --------------------- Technical Settings --------------------
		private int numIterations; // number of iterations needed to compute times
		// ------------------- Calc Method Parameters --------------------
		private Dictionary<int?, double[]> methodParams;

		/*
		 * this.methodParams[methodNum] = new Array(fa, ms, mv, is, iv);
		 *
		 * fa : fajr angle ms : maghrib selector (0 = angle; 1 = minutes after
		 * sunset) mv : maghrib parameter value (in angle or minutes) is : isha
		 * selector (0 = angle; 1 = minutes after maghrib) iv : isha parameter value
		 * (in angle or minutes)
		 */
		private double[] prayerTimesCurrent;
		private int[] offsets;

		public PrayTime()
		{
			// Initialize vars

			this.CalcMethod = 0;
			this.AsrJuristic = 0;
			this.DhuhrMinutes = 0;
			this.AdjustHighLats = 1;
			this.TimeFormat = 0;

			// Calculation Methods
			this.Jafari = 0; // Ithna Ashari
			this.Karachi = 1; // University of Islamic Sciences, Karachi
			this.ISNA = 2; // Islamic Society of North America (ISNA)
			this.MWL = 3; // Muslim World League (MWL)
			this.Makkah = 4; // Umm al-Qura, Makkah
			this.Egypt = 5; // Egyptian General Authority of Survey
			this.Tehran = 6; // Institute of Geophysics, University of Tehran
			this.Custom = 7; // Custom Setting

			// Juristic Methods
			this.Shafii = 0; // Shafii (standard)
			this.Hanafi = 1; // Hanafi

			// Adjusting Methods for Higher Latitudes
			this.None = 0; // No adjustment
			this.MidNight = 1; // middle of night
			this.OneSeventh = 2; // 1/7th of night
			this.AngleBased = 3; // angle/60th of night

			// Time Formats
			this.Time24 = 0; // 24-hour format
			this.Time12 = 1; // 12-hour format
			this.Time12NS = 2; // 12-hour format with no suffix
			this.Floating = 3; // floating point number

			// Time Names
			timeNames = new List<string>();
			timeNames.Add("Fajr");
			timeNames.Add("Sunrise");
			timeNames.Add("Dhuhr");
			timeNames.Add("Asr");
			timeNames.Add("Sunset");
			timeNames.Add("Maghrib");
			timeNames.Add("Isha");

			InvalidTime = "-----"; // The string used for invalid times

			// --------------------- Technical Settings --------------------

			this.NumIterations = 1; // number of iterations needed to compute
			// times

			// ------------------- Calc Method Parameters --------------------

			// Tuning offsets {fajr, sunrise, dhuhr, asr, sunset, maghrib, isha}
			offsets = new int[7];
			offsets[0] = 0;
			offsets[1] = 0;
			offsets[2] = 0;
			offsets[3] = 0;
			offsets[4] = 0;
			offsets[5] = 0;
			offsets[6] = 0;

			/*
			 *
			 * fa : fajr angle ms : maghrib selector (0 = angle; 1 = minutes after
			 * sunset) mv : maghrib parameter value (in angle or minutes) is : isha
			 * selector (0 = angle; 1 = minutes after maghrib) iv : isha parameter
			 * value (in angle or minutes)
			 */
			methodParams = new Dictionary<int?, double[]>();

			// Jafari
			double[] Jvalues = new double[] {16,0,4,0,14};
			methodParams[Convert.ToInt32(this.Jafari)] = Jvalues;

			// Karachi
			double[] Kvalues = new double[] {18,1,0,0,18};
			methodParams[Convert.ToInt32(this.Karachi)] = Kvalues;

			// ISNA
			double[] Ivalues = new double[] {15,1,0,0,15};
			methodParams[Convert.ToInt32(this.ISNA)] = Ivalues;

			// MWL
			double[] MWvalues = new double[] {18,1,0,0,17};
			methodParams[Convert.ToInt32(this.MWL)] = MWvalues;

			// Makkah
			double[] MKvalues = new double[] {18.5,1,0,1,90};
			methodParams[Convert.ToInt32(this.Makkah)] = MKvalues;

			// Egypt
			double[] Evalues = new double[] {19.5,1,0,0,17.5};
			methodParams[Convert.ToInt32(this.Egypt)] = Evalues;

			// Tehran
			double[] Tvalues = new double[] {17.7,0,4.5,0,14};
			methodParams[Convert.ToInt32(this.Tehran)] = Tvalues;

			// Custom
			double[] Cvalues = new double[] {18,1,0,0,17};
			methodParams[Convert.ToInt32(this.Custom)] = Cvalues;

		}

		// ---------------------- Trigonometric Functions -----------------------
		// range reduce angle in degrees.
		private double fixangle(double a)
		{

			a = a - (360 * (Math.Floor(a / 360.0)));

			a = a < 0 ? (a + 360) : a;

			return a;
		}

		// range reduce hours to 0..23
		private double fixhour(double a)
		{
			a = a - 24.0 * Math.Floor(a / 24.0);
			a = a < 0 ? (a + 24) : a;
			return a;
		}

		// radian to degree
		private double radiansToDegrees(double alpha)
		{
			return ((alpha * 180.0) / Math.PI);
		}

		// deree to radian
		private double DegreesToRadians(double alpha)
		{
			return ((alpha * Math.PI) / 180.0);
		}

		// degree sin
		private double dsin(double d)
		{
			return (Math.Sin(DegreesToRadians(d)));
		}

		// degree cos
		private double dcos(double d)
		{
			return (Math.Cos(DegreesToRadians(d)));
		}

		// degree tan
		private double dtan(double d)
		{
			return (Math.Tan(DegreesToRadians(d)));
		}

		// degree arcsin
		private double darcsin(double x)
		{
			double val = Math.Asin(x);
			return radiansToDegrees(val);
		}

		// degree arccos
		private double darccos(double x)
		{
			double val = Math.Acos(x);
			return radiansToDegrees(val);
		}

		// degree arctan
		private double darctan(double x)
		{
			double val = Math.Atan(x);
			return radiansToDegrees(val);
		}

		// degree arctan2
		private double darctan2(double y, double x)
		{
			double val = Math.Atan2(y, x);
			return radiansToDegrees(val);
		}

		// degree arccot
		private double darccot(double x)
		{
			double val = Math.Atan2(1.0, x);
			return radiansToDegrees(val);
		}

		// ---------------------- Time-Zone Functions -----------------------
		// compute local time-zone for a specific date
		private double TimeZone1
		{
			get
			{

               // TimeZone.CurrentTimeZone.IsDaylightSavingTime(new DateTime(year, month, day));

			//	TimeZone timez = TimeZone.Default;
			//	double hoursDiff = (timez.RawOffset / 1000.0) / 3600;
                return 1;// hoursDiff;
			}
		}

		// compute base time-zone of the system
		private double BaseTimeZone
		{
			get
			{
				//TimeZone timez = TimeZone.Default;
				//double hoursDiff = (timez.GetUtcOffset / 1000.0) / 3600;
                return 1;// hoursDiff;
    
			}
		}

		// detect daylight saving in a given date
		private double detectDaylightSaving()
		{
			//TimeZone timez = TimeZone.Default;
			//double hoursDiff = timez.DSTSavings;
            return 1;// hoursDiff;
		}

		// ---------------------- Julian Date Functions -----------------------
		// calculate julian date from a calendar date
		private double julianDate(int year, int month, int day)
		{

			if (month <= 2)
			{
				year -= 1;
				month += 12;
			}
			double A = Math.Floor(year / 100.0);

			double B = 2 - A + Math.Floor(A / 4.0);

			double JD = Math.Floor(365.25 * (year + 4716)) + Math.Floor(30.6001 * (month + 1)) + day + B - 1524.5;

			return JD;
		}

		// convert a calendar date to julian date (second method)
		private double calcJD(int year, int month, int day)
		{
			double J1970 = 2440588.0;
			DateTime date = new DateTime(year, month - 1, day);

			double ms = date.Ticks; // # of milliseconds since midnight Jan 1,
			// 1970
			double days = Math.Floor(ms / (1000.0 * 60.0 * 60.0 * 24.0));
			return J1970 + days - 0.5;

		}

		// ---------------------- Calculation Functions -----------------------
		// References:
		// http://www.ummah.net/astronomy/saltime
		// http://aa.usno.navy.mil/faq/docs/SunApprox.html
		// compute declination angle of sun and equation of time
		private double[] sunPosition(double jd)
		{

			double D = jd - 2451545;
			double g = fixangle(357.529 + 0.98560028 * D);
			double q = fixangle(280.459 + 0.98564736 * D);
			double L = fixangle(q + (1.915 * dsin(g)) + (0.020 * dsin(2 * g)));

			// double R = 1.00014 - 0.01671 * [self dcos:g] - 0.00014 * [self dcos:
			// (2*g)];
			double e = 23.439 - (0.00000036 * D);
			double d = darcsin(dsin(e) * dsin(L));
			double RA = (darctan2((dcos(e) * dsin(L)), (dcos(L)))) / 15.0;
			RA = fixhour(RA);
			double EqT = q / 15.0 - RA;
			double[] sPosition = new double[2];
			sPosition[0] = d;
			sPosition[1] = EqT;

			return sPosition;
		}

		// compute equation of time
		private double equationOfTime(double jd)
		{
			double eq = sunPosition(jd)[1];
			return eq;
		}

		// compute declination angle of sun
		private double sunDeclination(double jd)
		{
			double d = sunPosition(jd)[0];
			return d;
		}

		// compute mid-day (Dhuhr, Zawal) time
		private double computeMidDay(double t)
		{
			double T = equationOfTime(this.JDate + t);
			double Z = fixhour(12 - T);
			return Z;
		}

		// compute time for a given angle G
		private double computeTime(double G, double t)
		{

			double D = sunDeclination(this.JDate + t);
			double Z = computeMidDay(t);
			double Beg = -dsin(G) - dsin(D) * dsin(this.Lat);
			double Mid = dcos(D) * dcos(this.Lat);
			double V = darccos(Beg / Mid) / 15.0;

			return Z + (G > 90 ? - V : V);
		}

		// compute the time of Asr
		// Shafii: step=1, Hanafi: step=2
		private double computeAsr(double step, double t)
		{
			double D = sunDeclination(this.JDate + t);
			double G = -darccot(step + dtan(Math.Abs(this.Lat - D)));
			return computeTime(G, t);
		}

		// ---------------------- Misc Functions -----------------------
		// compute the difference between two times
		private double timeDiff(double time1, double time2)
		{
			return fixhour(time2 - time1);
		}

		// -------------------- Interface Functions --------------------
		// return prayer times for a given date
		private List<string> getDatePrayerTimes(int year, int month, int day, double latitude, double longitude, double tZone)
		{
			this.Lat = latitude;
			this.Lng = longitude;
			this.TimeZone = tZone;
			this.JDate = julianDate(year, month, day);
			double lonDiff = longitude / (15.0 * 24.0);
			this.JDate = this.JDate - lonDiff;
			return computeDayTimes();
		}

		// return prayer times for a given date
		private List<string> getPrayerTimes(DateTime date, double latitude, double longitude, double tZone)
		{

			int year = date.Year;
			int month = date.Month;
			int day = date.Day;

			return getDatePrayerTimes(year, month , day, latitude, longitude, tZone);
		}

		// set custom values for calculation parameters
		private double[] CustomParams
		{
			set
			{
    
				for (int i = 0; i < 5; i++)
				{
					if (value[i] == -1)
					{
						value[i] = methodParams[this.CalcMethod][i];
						methodParams[this.Custom] = value;
					}
					else
					{
						methodParams[this.Custom][i] = value[i];
					}
				}
				this.CalcMethod = this.Custom;
			}
		}

		// set the angle for calculating Fajr
		public virtual double FajrAngle
		{
			set
			{
				double[] @params = new double[] {value, -1, -1, -1, -1};
				CustomParams = @params;
			}
		}

		// set the angle for calculating Maghrib
		public virtual double MaghribAngle
		{
			set
			{
				double[] @params = new double[] {-1, 0, value, -1, -1};
				CustomParams = @params;
    
			}
		}

		// set the angle for calculating Isha
		public virtual double IshaAngle
		{
			set
			{
				double[] @params = new double[] {-1, -1, -1, 0, value};
				CustomParams = @params;
    
			}
		}

		// set the minutes after Sunset for calculating Maghrib
		public virtual double MaghribMinutes
		{
			set
			{
				double[] @params = new double[] {-1, 1, value, -1, -1};
				CustomParams = @params;
    
			}
		}

		// set the minutes after Maghrib for calculating Isha
		public virtual double IshaMinutes
		{
			set
			{
				double[] @params = new double[] {-1, -1, -1, 1, value};
				CustomParams = @params;
    
			}
		}

		// convert double hours to 24h format
		public virtual string floatToTime24(double time)
		{

			string result;

			if (double.IsNaN(time))
			{
				return InvalidTime;
			}

			time = fixhour(time + 0.5 / 60.0); // add 0.5 minutes to round
			int hours = (int)Math.Floor(time);
			double minutes = Math.Floor((time - hours) * 60.0);

			if ((hours >= 0 && hours <= 9) && (minutes >= 0 && minutes <= 9))
			{
				result = "0" + hours + ":0" + Math.Round(minutes);
			}
			else if ((hours >= 0 && hours <= 9))
			{
				result = "0" + hours + ":" + Math.Round(minutes);
			}
			else if ((minutes >= 0 && minutes <= 9))
			{
				result = hours + ":0" + Math.Round(minutes);
			}
			else
			{
				result = hours + ":" + Math.Round(minutes);
			}
			return result;
		}

		// convert double hours to 12h format
		public virtual string floatToTime12(double time, bool noSuffix)
		{

			if (double.IsNaN(time))
			{
				return InvalidTime;
			}

			time = fixhour(time + 0.5 / 60); // add 0.5 minutes to round
			int hours = (int)Math.Floor(time);
			double minutes = Math.Floor((time - hours) * 60);
			string suffix, result;
			if (hours >= 12)
			{
				suffix = "pm";
			}
			else
			{
				suffix = "am";
			}
			hours = ((((hours + 12) - 1) % (12)) + 1);
			/*hours = (hours + 12) - 1;
			int hrs = (int) hours % 12;
			hrs += 1;*/
			if (noSuffix == false)
			{
				if ((hours >= 0 && hours <= 9) && (minutes >= 0 && minutes <= 9))
				{
					result = "0" + hours + ":0" + Math.Round(minutes) + " " + suffix;
				}
				else if ((hours >= 0 && hours <= 9))
				{
					result = "0" + hours + ":" + Math.Round(minutes) + " " + suffix;
				}
				else if ((minutes >= 0 && minutes <= 9))
				{
					result = hours + ":0" + Math.Round(minutes) + " " + suffix;
				}
				else
				{
					result = hours + ":" + Math.Round(minutes) + " " + suffix;
				}

			}
			else
			{
				if ((hours >= 0 && hours <= 9) && (minutes >= 0 && minutes <= 9))
				{
					result = "0" + hours + ":0" + Math.Round(minutes);
				}
				else if ((hours >= 0 && hours <= 9))
				{
					result = "0" + hours + ":" + Math.Round(minutes);
				}
				else if ((minutes >= 0 && minutes <= 9))
				{
					result = hours + ":0" + Math.Round(minutes);
				}
				else
				{
					result = hours + ":" + Math.Round(minutes);
				}
			}
			return result;

		}

		// convert double hours to 12h format with no suffix
		public virtual string floatToTime12NS(double time)
		{
			return floatToTime12(time, true);
		}

		// ---------------------- Compute Prayer Times -----------------------
		// compute prayer times at given julian date
		private double[] computeTimes(double[] times)
		{

			double[] t = dayPortion(times);

			double Fajr = this.computeTime(180 - methodParams[this.CalcMethod][0], t[0]);

			double Sunrise = this.computeTime(180 - 0.833, t[1]);

			double Dhuhr = this.computeMidDay(t[2]);
			double Asr = this.computeAsr(1 + this.AsrJuristic, t[3]);
			double Sunset = this.computeTime(0.833, t[4]);

			double Maghrib = this.computeTime(methodParams[this.CalcMethod][2], t[5]);
			double Isha = this.computeTime(methodParams[this.CalcMethod][4], t[6]);

			double[] CTimes = new double[] {Fajr, Sunrise, Dhuhr, Asr, Sunset, Maghrib, Isha};

			return CTimes;

		}

		// compute prayer times at given julian date
		private List<string> computeDayTimes()
		{
			double[] times = new double[] {5, 6, 12, 13, 18, 18, 18}; // default times

			for (int i = 1; i <= this.NumIterations; i++)
			{
				times = computeTimes(times);
			}

			times = adjustTimes(times);
			times = tuneTimes(times);

			return adjustTimesFormat(times);
		}

		// adjust times in a prayer time array
		private double[] adjustTimes(double[] times)
		{
			for (int i = 0; i < times.Length; i++)
			{
				times[i] += this.TimeZone - this.Lng / 15;
			}

			times[2] += this.DhuhrMinutes / 60; // Dhuhr
			if (methodParams[this.CalcMethod][1] == 1) // Maghrib
			{
				times[5] = times[4] + methodParams[this.CalcMethod][2] / 60;
			}
			if (methodParams[this.CalcMethod][3] == 1) // Isha
			{
				times[6] = times[5] + methodParams[this.CalcMethod][4] / 60;
			}

			if (this.AdjustHighLats != this.None)
			{
				times = adjustHighLatTimes(times);
			}

			return times;
		}

		// convert times array to given time format
		private List<string> adjustTimesFormat(double[] times)
		{

			List<string> result = new List<string>();

			if (this.TimeFormat == this.Floating)
			{
				foreach (double time in times)
				{
					result.Add(Convert.ToString(time));
				}
				return result;
			}

			for (int i = 0; i < 7; i++)
			{
				if (this.TimeFormat == this.Time12)
				{
					result.Add(floatToTime12(times[i], false));
				}
				else if (this.TimeFormat == this.Time12NS)
				{
					result.Add(floatToTime12(times[i], true));
				}
				else
				{
					result.Add(floatToTime24(times[i]));
				}
			}
			return result;
		}

		// adjust Fajr, Isha and Maghrib for locations in higher latitudes
		private double[] adjustHighLatTimes(double[] times)
		{
			double nightTime = timeDiff(times[4], times[1]); // sunset to sunrise

			// Adjust Fajr
			double FajrDiff = nightPortion(methodParams[this.CalcMethod][0]) * nightTime;

			if (double.IsNaN(times[0]) || timeDiff(times[0], times[1]) > FajrDiff)
			{
				times[0] = times[1] - FajrDiff;
			}

			// Adjust Isha
			double IshaAngle = (methodParams[this.CalcMethod][3] == 0) ? methodParams[this.CalcMethod][4] : 18;
			double IshaDiff = this.nightPortion(IshaAngle) * nightTime;
			if (double.IsNaN(times[6]) || this.timeDiff(times[4], times[6]) > IshaDiff)
			{
				times[6] = times[4] + IshaDiff;
			}

			// Adjust Maghrib
			double MaghribAngle = (methodParams[this.CalcMethod][1] == 0) ? methodParams[this.CalcMethod][2] : 4;
			double MaghribDiff = nightPortion(MaghribAngle) * nightTime;
			if (double.IsNaN(times[5]) || this.timeDiff(times[4], times[5]) > MaghribDiff)
			{
				times[5] = times[4] + MaghribDiff;
			}

			return times;
		}

		// the night portion used for adjusting times in higher latitudes
		private double nightPortion(double angle)
		{
		   double calc = 0;

		if (adjustHighLats == AngleBased_Renamed)
		{
			calc = (angle) / 60.0;
		}
		else if (adjustHighLats == MidNight_Renamed)
		{
			calc = 0.5;
		}
		else if (adjustHighLats == OneSeventh_Renamed)
		{
			calc = 0.14286;
		}

		return calc;
		}

		// convert hours to day portions
		private double[] dayPortion(double[] times)
		{
			for (int i = 0; i < 7; i++)
			{
				times[i] /= 24;
			}
			return times;
		}

		// Tune timings for adjustments
		// Set time offsets
		public virtual void tune(int[] offsetTimes)
		{

			for (int i = 0; i < offsetTimes.Length; i++) // offsetTimes length
			{
				// should be 7 in order
				// of Fajr, Sunrise,
				// Dhuhr, Asr, Sunset,
				// Maghrib, Isha
				this.offsets[i] = offsetTimes[i];
			}
		}

		private double[] tuneTimes(double[] times)
		{
			for (int i = 0; i < times.Length; i++)
			{
				times[i] = times[i] + this.offsets[i] / 60.0;
			}

			return times;
		}

		/// <param name="args"> </param>
		public static void Main(string[] args)
		{
			double latitude = 24.8351823817939;
			double longitude = 46.758465366438;
			double timezone = 3;
			// Test Prayer times here
			PrayTime prayers = new PrayTime();

			prayers.TimeFormat = prayers.Time12_Renamed;
			prayers.CalcMethod = prayers.Makkah_Renamed;
			prayers.AsrJuristic = prayers.Shafii_Renamed;
			prayers.AdjustHighLats = 1;
			int[] offsets = new int[] {0, 0, 0, 0, 0, 0, 0}; // {Fajr,Sunrise,Dhuhr,Asr,Sunset,Maghrib,Isha}
			prayers.tune(offsets);

			DateTime now = DateTime.Now;
			DateTime cal = new DateTime();
			cal = now;

			List<string> prayerTimes = prayers.getPrayerTimes(cal, latitude, longitude, timezone);
			List<string> prayerNames = prayers.TimeNames;

			for (int i = 0; i < prayerTimes.Count; i++)
			{
				Console.WriteLine(prayerNames[i] + " - " + prayerTimes[i]);
			}

		}

		public virtual int CalcMethod
		{
			get
			{
				return calcMethod;
			}
			set
			{
				this.calcMethod = value;
			}
		}


		public virtual int AsrJuristic
		{
			get
			{
				return asrJuristic;
			}
			set
			{
				this.asrJuristic = value;
			}
		}


		public virtual int DhuhrMinutes
		{
			get
			{
				return dhuhrMinutes;
			}
			set
			{
				this.dhuhrMinutes = value;
			}
		}


		public virtual int AdjustHighLats
		{
			get
			{
				return adjustHighLats;
			}
			set
			{
				this.adjustHighLats = value;
			}
		}


		public virtual int TimeFormat
		{
			get
			{
				return timeFormat;
			}
			set
			{
				this.timeFormat = value;
			}
		}


		public virtual double Lat
		{
			get
			{
				return lat;
			}
			set
			{
				this.lat = value;
			}
		}


		public virtual double Lng
		{
			get
			{
				return lng;
			}
			set
			{
				this.lng = value;
			}
		}


		public virtual double TimeZone
		{
			get
			{
				return timeZone;
			}
			set
			{
				this.timeZone = value;
			}
		}


		public virtual double JDate
		{
			get
			{
				return JDate_Renamed;
			}
			set
			{
				JDate_Renamed = value;
			}
		}


		private int Jafari
		{
			get
			{
				return Jafari_Renamed;
			}
			set
			{
				Jafari_Renamed = value;
			}
		}


		private int Karachi
		{
			get
			{
				return Karachi_Renamed;
			}
			set
			{
				Karachi_Renamed = value;
			}
		}


		private int ISNA
		{
			get
			{
				return ISNA_Renamed;
			}
			set
			{
				ISNA_Renamed = value;
			}
		}


		private int MWL
		{
			get
			{
				return MWL_Renamed;
			}
			set
			{
				MWL_Renamed = value;
			}
		}


		private int Makkah
		{
			get
			{
				return Makkah_Renamed;
			}
			set
			{
				Makkah_Renamed = value;
			}
		}


		private int Egypt
		{
			get
			{
				return Egypt_Renamed;
			}
			set
			{
				Egypt_Renamed = value;
			}
		}


		private int Custom
		{
			get
			{
				return Custom_Renamed;
			}
			set
			{
				Custom_Renamed = value;
			}
		}


		private int Tehran
		{
			get
			{
				return Tehran_Renamed;
			}
			set
			{
				Tehran_Renamed = value;
			}
		}


		private int Shafii
		{
			get
			{
				return Shafii_Renamed;
			}
			set
			{
				Shafii_Renamed = value;
			}
		}


		private int Hanafi
		{
			get
			{
				return Hanafi_Renamed;
			}
			set
			{
				Hanafi_Renamed = value;
			}
		}


		private int None
		{
			get
			{
				return None_Renamed;
			}
			set
			{
				None_Renamed = value;
			}
		}


		private int MidNight
		{
			get
			{
				return MidNight_Renamed;
			}
			set
			{
				MidNight_Renamed = value;
			}
		}


		private int OneSeventh
		{
			get
			{
				return OneSeventh_Renamed;
			}
			set
			{
				OneSeventh_Renamed = value;
			}
		}


		private int AngleBased
		{
			get
			{
				return AngleBased_Renamed;
			}
			set
			{
				AngleBased_Renamed = value;
			}
		}


		private int Time24
		{
			get
			{
				return Time24_Renamed;
			}
			set
			{
				Time24_Renamed = value;
			}
		}


		private int Time12
		{
			get
			{
				return Time12_Renamed;
			}
			set
			{
				Time12_Renamed = value;
			}
		}


		private int Time12NS
		{
			get
			{
				return Time12NS_Renamed;
			}
			set
			{
				Time12NS_Renamed = value;
			}
		}


		private int Floating
		{
			get
			{
				return Floating_Renamed;
			}
			set
			{
				Floating_Renamed = value;
			}
		}


		private int NumIterations
		{
			get
			{
				return numIterations;
			}
			set
			{
				this.numIterations = value;
			}
		}


		public virtual List<string> TimeNames
		{
			get
			{
				return timeNames;
			}
		}
	}

}