//   libPolinpol.h
//
//   Library for polynomial interpolation, using Catmull-Rom and/or the
//   Three-Point-Spline, both uniform and recursive version.
//   See also my paper "Cubic Spline Interpolation in Real-Time Applications
//   using Three Control Points" published in the Proceedings of the 27.
//   International Conference in Central Europe on Computer Graphics, Visualization
//   and Computer Vision'2019, vol 2, pp 1-10, online version:
//   http://wscg.zcu.cz/WSCG2019/!!_CSRN-2901.pdf
//
//   If you need performance, you might want to calculate the s-values in advance.
//   Omitted here for simplicity.
//   I use floats since it is aimed at embedded systems and computer graphics, but
//   you can replace them without any problems for doubles if you need the
//   additional precision.
//   Talking of which: the way I define the recursive functions (i.e. the s-values
//   based on one segment) is better precisionwise than the more common definition
//   (i.e. the s-values based on a sum of segments) - the latter will loose precision
//   after many segments, especially if there are many long segments followed by a
//   short one.
//
//   A word on timings, which is important to keep the mathematical properties.
//   For uniform interpolation, each segment needs to take the exact same length
//   of time.
//   For recursive, the time for each segment is dependent on the s-value
//   representing it, e.g. f*sr for segment r, with f a scalefactor that needs to
//   be the same for each segment, during the whole interpolation.
//
//   If you have any further questions/comments, feel free to contact me under
//   j.ogniewski@gmail.com. I will do my best to answer as soon as I can.
//
//   Copyright 2019 Jens Ogniewski
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the full License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//                                 Apache License
//                           Version 2.0, January 2004
//                        http://www.apache.org/licenses/
//
//   TERMS AND CONDITIONS FOR USE, REPRODUCTION, AND DISTRIBUTION
//
//   1. Definitions.
//
//      "License" shall mean the terms and conditions for use, reproduction,
//      and distribution as defined by Sections 1 through 9 of this document.
//
//      "Licensor" shall mean the copyright owner or entity authorized by
//      the copyright owner that is granting the License.
//
//      "Legal Entity" shall mean the union of the acting entity and all
//      other entities that control, are controlled by, or are under common
//      control with that entity. For the purposes of this definition,
//      "control" means (i) the power, direct or indirect, to cause the
//      direction or management of such entity, whether by contract or
//      otherwise, or (ii) ownership of fifty percent (50%) or more of the
//      outstanding shares, or (iii) beneficial ownership of such entity.
//
//      "You" (or "Your") shall mean an individual or Legal Entity
//      exercising permissions granted by this License.
//
//      "Source" form shall mean the preferred form for making modifications,
//      including but not limited to software source code, documentation
//      source, and configuration files.

//      "Object" form shall mean any form resulting from mechanical
//      transformation or translation of a Source form, including but
//      not limited to compiled object code, generated documentation,
//      and conversions to other media types.
//
//      "Work" shall mean the work of authorship, whether in Source or
//      Object form, made available under the License, as indicated by a
//      copyright notice that is included in or attached to the work
//      (an example is provided in the Appendix below).
//
//      "Derivative Works" shall mean any work, whether in Source or Object
//      form, that is based on (or derived from) the Work and for which the
//      editorial revisions, annotations, elaborations, or other modifications
//      represent, as a whole, an original work of authorship. For the purposes
//      of this License, Derivative Works shall not include works that remain
//      separable from, or merely link (or bind by name) to the interfaces of,
//      the Work and Derivative Works thereof.
//
//      "Contribution" shall mean any work of authorship, including
//      the original version of the Work and any modifications or additions
//      to that Work or Derivative Works thereof, that is intentionally
//      submitted to Licensor for inclusion in the Work by the copyright owner
//      or by an individual or Legal Entity authorized to submit on behalf of
//      the copyright owner. For the purposes of this definition, "submitted"
//      means any form of electronic, verbal, or written communication sent
//      to the Licensor or its representatives, including but not limited to
//      communication on electronic mailing lists, source code control systems,
//      and issue tracking systems that are managed by, or on behalf of, the
//      Licensor for the purpose of discussing and improving the Work, but
//      excluding communication that is conspicuously marked or otherwise
//      designated in writing by the copyright owner as "Not a Contribution."
//
//      "Contributor" shall mean Licensor and any individual or Legal Entity
//      on behalf of whom a Contribution has been received by Licensor and
//      subsequently incorporated within the Work.
//
//   2. Grant of Copyright License. Subject to the terms and conditions of
//      this License, each Contributor hereby grants to You a perpetual,
//      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
//      copyright license to reproduce, prepare Derivative Works of,
//      publicly display, publicly perform, sublicense, and distribute the
//      Work and such Derivative Works in Source or Object form.
//
//   3. Grant of Patent License. Subject to the terms and conditions of
//      this License, each Contributor hereby grants to You a perpetual,
//      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
//      (except as stated in this section) patent license to make, have made,
//      use, offer to sell, sell, import, and otherwise transfer the Work,
//      where such license applies only to those patent claims licensable
//      by such Contributor that are necessarily infringed by their
//      Contribution(s) alone or by combination of their Contribution(s)
//      with the Work to which such Contribution(s) was submitted. If You
//      institute patent litigation against any entity (including a
//      cross-claim or counterclaim in a lawsuit) alleging that the Work
//      or a Contribution incorporated within the Work constitutes direct
//      or contributory patent infringement, then any patent licenses
//      granted to You under this License for that Work shall terminate
//      as of the date such litigation is filed.
//
//   4. Redistribution. You may reproduce and distribute copies of the
//      Work or Derivative Works thereof in any medium, with or without
//      modifications, and in Source or Object form, provided that You
//      meet the following conditions:
//
//      (a) You must give any other recipients of the Work or
//          Derivative Works a copy of this License; and
//
//      (b) You must cause any modified files to carry prominent notices
//          stating that You changed the files; and
//
//      (c) You must retain, in the Source form of any Derivative Works
//          that You distribute, all copyright, patent, trademark, and
//          attribution notices from the Source form of the Work,
//          excluding those notices that do not pertain to any part of
//          the Derivative Works; and
//
//      (d) If the Work includes a "NOTICE" text file as part of its
//          distribution, then any Derivative Works that You distribute must
//          include a readable copy of the attribution notices contained
//          within such NOTICE file, excluding those notices that do not
//          pertain to any part of the Derivative Works, in at least one
//          of the following places: within a NOTICE text file distributed
//          as part of the Derivative Works; within the Source form or
//          documentation, if provided along with the Derivative Works; or,
//          within a display generated by the Derivative Works, if and
//          wherever such third-party notices normally appear. The contents
//          of the NOTICE file are for informational purposes only and
//          do not modify the License. You may add Your own attribution
//          notices within Derivative Works that You distribute, alongside
//          or as an addendum to the NOTICE text from the Work, provided
//          that such additional attribution notices cannot be construed
//          as modifying the License.
//
//      You may add Your own copyright statement to Your modifications and
//      may provide additional or different license terms and conditions
//      for use, reproduction, or distribution of Your modifications, or
//      for any such Derivative Works as a whole, provided Your use,
//      reproduction, and distribution of the Work otherwise complies with
//      the conditions stated in this License.
//
//   5. Submission of Contributions. Unless You explicitly state otherwise,
//      any Contribution intentionally submitted for inclusion in the Work
//      by You to the Licensor shall be under the terms and conditions of
//      this License, without any additional terms or conditions.
//      Notwithstanding the above, nothing herein shall supersede or modify
//      the terms of any separate license agreement you may have executed
//      with Licensor regarding such Contributions.
//
//   6. Trademarks. This License does not grant permission to use the trade
//      names, trademarks, service marks, or product names of the Licensor,
//      except as required for reasonable and customary use in describing the
//      origin of the Work and reproducing the content of the NOTICE file.
//
//   7. Disclaimer of Warranty. Unless required by applicable law or
//      agreed to in writing, Licensor provides the Work (and each
//      Contributor provides its Contributions) on an "AS IS" BASIS,
//      WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
//      implied, including, without limitation, any warranties or conditions
//      of TITLE, NON-INFRINGEMENT, MERCHANTABILITY, or FITNESS FOR A
//      PARTICULAR PURPOSE. You are solely responsible for determining the
//      appropriateness of using or redistributing the Work and assume any
//      risks associated with Your exercise of permissions under this License.
//
//   8. Limitation of Liability. In no event and under no legal theory,
//      whether in tort (including negligence), contract, or otherwise,
//      unless required by applicable law (such as deliberate and grossly
//      negligent acts) or agreed to in writing, shall any Contributor be
//      liable to You for damages, including any direct, indirect, special,
//      incidental, or consequential damages of any character arising as a
//      result of this License or out of the use or inability to use the
//      Work (including but not limited to damages for loss of goodwill,
//      work stoppage, computer failure or malfunction, or any and all
//      other commercial damages or losses), even if such Contributor
//      has been advised of the possibility of such damages.
//
//   9. Accepting Warranty or Additional Liability. While redistributing
//      the Work or Derivative Works thereof, You may choose to offer,
//      and charge a fee for, acceptance of support, warranty, indemnity,
//      or other liability obligations and/or rights consistent with this
//      License. However, in accepting such obligations, You may act only
//      on Your own behalf and on Your sole responsibility, not on behalf
//      of any other Contributor, and only if You agree to indemnify,
//      defend, and hold each Contributor harmless for any liability
//      incurred by, or claims asserted against, such Contributor by reason
//      of your accepting any such warranty or additional liability.
//
//   END OF TERMS AND CONDITIONS
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.
//


void uniformCrInterpolation(float u, float alpha, float *controlPoints, int segment, float *respos, int dimensions)
{
	int i;
	float ipd0=alpha*u*(u-1.0)*(u-1.0);
	float ipf =u*u*(3.0-2.0*u);
	float ipd1=alpha*u*u*(u-1.0);

	for (i=0; i<dimensions; ++i) {
		respos[i]=controlPoints[(segment)*dimensions+i];
		if (controlPoints[(segment+1)*dimensions+i]!=controlPoints[(segment-1)*dimensions+i]) {
			respos[i]=respos[i]+ipd0*(controlPoints[(segment+1)*dimensions+i]-controlPoints[(segment-1)*dimensions+i]);
		}
		respos[i]=respos[i]   +ipf *(controlPoints[(segment+1)*dimensions+i]-controlPoints[(segment)*dimensions+i]);
		if (controlPoints[(segment+2)*dimensions+i]!=controlPoints[(segment)*dimensions+i]) {
			respos[i]=respos[i]+ipd1*(controlPoints[(segment+2)*dimensions+i]-controlPoints[(segment)*dimensions+i]);
		}
	}
}

void uniformTpsInterpolation(float u, float alpha, float *controlPoints, int segment, float *respos, int dimensions)
{
	int i;
	float ipd0=alpha*u*(u-1.0)*(u-1.0);
	float ipf =u*u*(3.0-2.0*u);
	float ipd1=alpha*u*u*(u-1.0);

	for (i=0; i<dimensions; ++i) {
		respos[i]=controlPoints[(segment)*dimensions+i];
		if (controlPoints[(segment)*dimensions+i]!=controlPoints[(segment-1)*dimensions+i]) {
			respos[i]=respos[i]+ipd0*(controlPoints[(segment)*dimensions+i]-controlPoints[(segment-1)*dimensions+i]);
		}
		respos[i]=respos[i]+ipf    *(controlPoints[(segment+1)*dimensions+i]-controlPoints[(segment)*dimensions+i]);
		respos[i]=respos[i]+ipd1   *(controlPoints[(segment+1)*dimensions+i]-controlPoints[(segment)*dimensions+i]);
	}
}

//note: segment = r;
void recursiveTpsInterpolation(float v, float beta, float *controlPoints, int segment, float *respos, int dimensions)
{
	int i;
	float srm1=0.0;
	float sr  =0.0;
	float v1=0.0;
	float v2=0.0;
	float v3=1.0;

	for (i=0; i<dimensions; ++i) {
		srm1=srm1+(controlPoints[(segment)*dimensions+i]-controlPoints[(segment-1)*dimensions+i])
		         *(controlPoints[(segment)*dimensions+i]-controlPoints[(segment-1)*dimensions+i]);
		sr  =sr  +(controlPoints[(segment+1)*dimensions+i]-controlPoints[(segment)*dimensions+i])
		         *(controlPoints[(segment+1)*dimensions+i]-controlPoints[(segment)*dimensions+i]);
	}
	srm1=sqrt(srm1);
	sr  =sqrt(sr);
	srm1=pow(srm1,beta);
	sr  =pow(sr,beta);

	if (sr>0.0) {
		v1=v/sr;
		v3=1.0-v1;
	}

	if (srm1>0.0) {
		v2=v/srm1;
	}
	for (i=0; i<dimensions; ++i) {
		respos[i]= controlPoints[segment*dimensions+i]
		          +v3*v3*v2*     (controlPoints[segment*dimensions+i]-controlPoints[(segment-1)*dimensions+i])
		          +v1*(v1+v1*v3)*(controlPoints[(segment+1)*dimensions+i]-controlPoints[segment*dimensions+i]);
	}
}

//note: segment = r;
void recursiveCrInterpolation(float v, float beta, float *controlPoints, int segment, float *respos, int dimensions)
{
	int i;
	float srm1=0.0;
	float sr  =0.0;
	float srp1=0.0;

	for (i=0; i<dimensions; ++i) {
		srm1=srm1+(controlPoints[(segment)*dimensions+i]-controlPoints[(segment-1)*dimensions+i])
		         *(controlPoints[(segment)*dimensions+i]-controlPoints[(segment-1)*dimensions+i]);
		sr  =sr  +(controlPoints[(segment+1)*dimensions+i]-controlPoints[(segment)*dimensions+i])
		         *(controlPoints[(segment+1)*dimensions+i]-controlPoints[(segment)*dimensions+i]);
		srp1=srp1+(controlPoints[(segment+2)*dimensions+i]-controlPoints[(segment+1)*dimensions+i])
		         *(controlPoints[(segment+2)*dimensions+i]-controlPoints[(segment+1)*dimensions+i]);
	}
	srm1=sqrt(srm1);
	sr  =sqrt(sr);
	srp1=sqrt(srp1);
	srm1=pow(srm1,beta);
	sr  =pow(sr,beta);
	srp1=pow(srp1,beta);

	float v1=0.0;
	float v2=0.0;
	float v3=(v-sr);
	float v4=v;
	float omv=1.0-v1;
	float tl01=srm1+sr;

	if (sr>0.0) {
		v1=v/sr;
		omv=1.0-v1;
	}

	if (srm1>0.0) {
		v2=v/srm1;
	}

	if (srp1>0.0) {
		v3=v3/srp1;
	}

	if (sr+srp1>0.0) {
		v4=v4/(sr+srp1);
	}

	for (i=0; i<dimensions; ++i) {
		float tp11,tp12,tp13,tp21,tp22;

		tp11=-v2*controlPoints[(segment-1)*dimensions+i]     +(1.0+v2)*controlPoints[segment*dimensions+i];
		tp12=omv*controlPoints[segment*dimensions+i]         +v1*controlPoints[(segment+1)*dimensions+i];
		tp13=(1.0-v3)*controlPoints[(segment+1)*dimensions+i]+v3*controlPoints[(segment+2)*dimensions+i];

		tp21=(sr-v)/tl01*tp11+(srm1+v)/tl01*tp12;
		tp22=(1.0-v4)*tp12+v4*tp13;
		
		respos[i]=omv*tp21+v1*tp22;
	}
}

//in case different betas should be used, e.g. for change of control-points
//note: segment = r;
void recursiveTpsInterpolationDifferentBetas(float v, float *betas, float *controlPoints, int segment, float *respos, int dimensions)
{
	int i;
	float srm1=0.0;
	float sr  =0.0;
	float v1=0.0;
	float v2=0.0;
	float v3=1.0;

	for (i=0; i<dimensions; ++i) {
		srm1=srm1+(controlPoints[(segment)*dimensions+i]-controlPoints[(segment-1)*dimensions+i])
		         *(controlPoints[(segment)*dimensions+i]-controlPoints[(segment-1)*dimensions+i]);
		sr  =sr  +(controlPoints[(segment+1)*dimensions+i]-controlPoints[(segment)*dimensions+i])
		         *(controlPoints[(segment+1)*dimensions+i]-controlPoints[(segment)*dimensions+i]);
	}
	srm1=sqrt(srm1);
	sr  =sqrt(sr);
	srm1=pow(srm1,betas[segment-1]);
	sr  =pow(sr,betas[segment]);

	if (sr>0.0) {
		v1=v/sr;
		v3=1.0-v1;
	}

	if (srm1>0.0) {
		v2=v/srm1;
	}
	for (i=0; i<dimensions; ++i) {
		respos[i]= controlPoints[segment*dimensions+i]
		          +v3*v3*v2*     (controlPoints[segment*dimensions+i]-controlPoints[(segment-1)*dimensions+i])
		          +v1*(v1+v1*v3)*(controlPoints[(segment+1)*dimensions+i]-controlPoints[segment*dimensions+i]);
	}
}

//note: segment = r;
void recursiveCrInterpolationDifferentBetas(float v, float *betas, float *controlPoints, int segment, float *respos, int dimensions)
{
	int i;
	float srm1=0.0;
	float sr  =0.0;
	float srp1=0.0;

	for (i=0; i<dimensions; ++i) {
		srm1=srm1+(controlPoints[(segment)*dimensions+i]-controlPoints[(segment-1)*dimensions+i])
		         *(controlPoints[(segment)*dimensions+i]-controlPoints[(segment-1)*dimensions+i]);
		sr  =sr  +(controlPoints[(segment+1)*dimensions+i]-controlPoints[(segment)*dimensions+i])
		         *(controlPoints[(segment+1)*dimensions+i]-controlPoints[(segment)*dimensions+i]);
		srp1=srp1+(controlPoints[(segment+2)*dimensions+i]-controlPoints[(segment+1)*dimensions+i])
		         *(controlPoints[(segment+2)*dimensions+i]-controlPoints[(segment+1)*dimensions+i]);
	}
	srm1=sqrt(srm1);
	sr  =sqrt(sr);
	srp1=sqrt(srp1);
	srm1=pow(srm1,betas[segment-1]);
	sr  =pow(sr,betas[segment]);
	srp1=pow(srp1,betas[segment+1]);

	float v1=0.0;
	float v2=0.0;
	float v3=(v-sr);
	float v4=v;
	float omv=1.0-v1;
	float tl01=srm1+sr;

	if (sr>0.0) {
		v1=v/sr;
		omv=1.0-v1;
	}

	if (srm1>0.0) {
		v2=v/srm1;
	}

	if (srp1>0.0) {
		v3=v3/srp1;
	}

	if (sr+srp1>0.0) {
		v4=v4/(sr+srp1);
	}

	for (i=0; i<dimensions; ++i) {
		float tp11,tp12,tp13,tp21,tp22;

		tp11=-v2*controlPoints[(segment-1)*dimensions+i]     +(1.0+v2)*controlPoints[segment*dimensions+i];
		tp12=omv*controlPoints[segment*dimensions+i]         +v1*controlPoints[(segment+1)*dimensions+i];
		tp13=(1.0-v3)*controlPoints[(segment+1)*dimensions+i]+v3*controlPoints[(segment+2)*dimensions+i];

		tp21=(sr-v)/tl01*tp11+(srm1+v)/tl01*tp12;
		tp22=(1.0-v4)*tp12+v4*tp13;
		
		respos[i]=omv*tp21+v1*tp22;
	}
}

void recursiveCrDerivativeInPoint(float v, float beta, float *controlPoints, int segment, float *respos, int dimensions) {
	int i;
	float srm1=0.0;
	float sr  =0.0;
	float srp1=0.0;

	for (i=0; i<dimensions; ++i) {
		srm1=srm1+(controlPoints[(segment)*dimensions+i]-controlPoints[(segment-1)*dimensions+i])
		         *(controlPoints[(segment)*dimensions+i]-controlPoints[(segment-1)*dimensions+i]);
		sr  =sr  +(controlPoints[(segment+1)*dimensions+i]-controlPoints[(segment)*dimensions+i])
		         *(controlPoints[(segment+1)*dimensions+i]-controlPoints[(segment)*dimensions+i]);
		srp1=srp1+(controlPoints[(segment+2)*dimensions+i]-controlPoints[(segment+1)*dimensions+i])
		         *(controlPoints[(segment+2)*dimensions+i]-controlPoints[(segment+1)*dimensions+i]);
	}
	srm1=sqrt(srm1);
	sr  =sqrt(sr);
	srp1=sqrt(srp1);
	srm1=pow(srm1,beta);
	sr  =pow(sr,beta);
	srp1=pow(srp1,beta);

	float s0Ps1 = srm1+sr;
	float s1Ms0Ps1 = sr*s0Ps1;
	float s1Ms1Ps2 = sr*(sr+srp1);

	for (i=0; i<dimensions; ++i) {
		float bw=0.0;
		if (sr*s1Ms1Ps2>0.0) {
			bw=(2.0*sr*v-3.0*v*v)/(sr*s1Ms1Ps2);
		}
		if (sr*s1Ms0Ps1>0.0) {
			bw=bw+(srm1*sr+4.0*sr*v-3.0*v*v)/(sr*s1Ms0Ps1);
		}
		respos[i] = bw*(controlPoints[(segment+1)*dimensions+i]-controlPoints[(segment)*dimensions+i]);
		if (srp1*s1Ms1Ps2>0.0) {
			respos[i] = respos[i] + (3.0*v*v-2.0*v*sr)/(srp1*s1Ms1Ps2)*(controlPoints[(segment+2)*dimensions+i]-controlPoints[(segment+1)*dimensions+i]);
		}
		if (srm1*s1Ms0Ps1>0.0) {
			respos[i] = respos[i] + (sr*sr-sr*4.0*v+3.0*v*v)/(srm1*s1Ms0Ps1)*(controlPoints[(segment)*dimensions+i]-controlPoints[(segment-1)*dimensions+i]);
		}
	}
}

void recursiveTpsDerivativeInPoint(float v, float beta, float *controlPoints, int segment, float *respos, int dimensions) {
	int i;
	float srm1=0.0;
	float sr  =0.0;
	float tv  =0.0;
	float bw  =0.0;

	for (i=0; i<dimensions; ++i) {
		srm1=srm1+(controlPoints[(segment)*dimensions+i]-controlPoints[(segment-1)*dimensions+i])
		         *(controlPoints[(segment)*dimensions+i]-controlPoints[(segment-1)*dimensions+i]);
		sr  =sr  +(controlPoints[(segment+1)*dimensions+i]-controlPoints[(segment)*dimensions+i])
		         *(controlPoints[(segment+1)*dimensions+i]-controlPoints[(segment)*dimensions+i]);
	}

	srm1=sqrt(srm1);
	sr  =sqrt(sr);
	srm1=pow(srm1,beta);
	sr  =pow(sr,beta);

	if (sr>0.0) {
		tv=v/sr;
		bw=4.0*tv-3.0*tv*tv;
	}

	for (i=0; i<dimensions; ++i) {
		respos[i]=0.0;
		if (sr>0.0) {
			respos[i]=respos[i]+bw/sr*(controlPoints[(segment+1)*dimensions+i]-controlPoints[segment*dimensions+i]);
		}
		if (srm1>0.0) {
			respos[i]=respos[i]+(1.0-bw)/srm1*(controlPoints[(segment)*dimensions+i]-controlPoints[(segment-1)*dimensions+i]);
		}
	}
}

//NOTE if you want to change interpolation points, you either want to a) insert a point or b) replace the next point (and in
//any case the current point to keep the continuity)
//Case a) is describe in the paper and in the demo, b) is used here for simplicity.
//Also, note that it will replace the current point, the next point and the last point (the latter is needed to keep c1 continuity).

void changeControlPointsRecursiveTps(float v, float beta, float *controlPoints, float *evasionPoint, int segment, int dimensions) {
	int i;
	float newPoint[dimensions];
	float derivativeInNewPoint[dimensions];

	recursiveTpsInterpolation(v, beta, controlPoints, segment, newPoint, dimensions);
	recursiveTpsDerivativeInPoint(v, beta, controlPoints, segment, derivativeInNewPoint, dimensions);

	for (i=0; i<dimensions; ++i) {
		controlPoints[(segment-1)*dimensions+i]=newPoint[i]-derivativeInNewPoint[i];
		controlPoints[(segment)*dimensions+i]=newPoint[i];
		controlPoints[(segment+1)*dimensions+i]=evasionPoint[i];
	}
}

void changeControlPointsRecursiveCr(float v, float beta, float *controlPoints, float *evasionPoint, int segment, int dimensions) {
	int i;
	float newsegment=0.0;
	float newPoint[dimensions];
	float derivativeInNewPoint[dimensions];

	recursiveCrInterpolation(v, beta, controlPoints, segment, newPoint, dimensions);
	recursiveCrDerivativeInPoint(v, beta, controlPoints, segment, derivativeInNewPoint, dimensions);

	for (i=0; i<dimensions; ++i) {
		newsegment=newsegment+(evasionPoint[(segment)*dimensions+i]-newPoint[i])
		                       *(evasionPoint[(segment)*dimensions+i]-newPoint[i]);
	}
	newsegment=sqrt(newsegment);
	newsegment=pow(newsegment,beta);

	for (i=0; i<dimensions; ++i) {
		controlPoints[(segment-1)*dimensions+i]=0.0;
		if (newsegment>0.0) {
			derivativeInNewPoint[i]= ( (1.0+newsegment)*derivativeInNewPoint[i]
			                          -((evasionPoint[(segment)*dimensions+i]-newPoint[(segment-1)*dimensions+i])/newsegment))
			                         /newsegment;
		}
		controlPoints[(segment-1)*dimensions+i]=newPoint[i]-derivativeInNewPoint[i];
		controlPoints[(segment)*dimensions+i]  =newPoint[i];
		controlPoints[(segment+1)*dimensions+i]=evasionPoint[i];
	}
}

