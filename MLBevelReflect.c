/*
	GIMP plugin for version 2.8
	Copyright (C) 2013, Marcel Lancelle,  gimp <at> marcel-lancelle.de

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#define PLUGIN_VERSION "2013-03-19"
#define PLUGIN_PROC "plug-in-mlbevelreflect"

#include <libgimp/gimp.h>
#include <libgimp/gimpui.h>

enum SlopeTypeEnum{SLOPE_FLAT=1, SLOPE_ROUND, SLOPE_ROUND2}; // 1: flat, 2: round and smooth
typedef struct {
	gboolean invertInput; // if false use white objects on black background
	gint slopeType;
	gint maxNormalMapBlurRadius; // for slope type 1
	gint probeImgDrawable; // id of drawable, -1 for normal map
	gdouble specularity; // 0: diffuse, 1: fully specular
	gboolean preview;
} MLBevReflVals;

struct { // copies of pointers to widgets, for access in callbacks
	GimpPreview* preview;
	GtkCheckButton* invertInput_checkButton;
	GtkRadioButton* slope_flat_radiobutton;
	GtkRadioButton* slope_round_radiobutton;
	GtkRadioButton* slope_round2_radiobutton;
	GtkRadioButton* probeImg_radiobutton;
	GtkRadioButton* normalMap_radiobutton;
	GimpDrawableComboBox* probeImg_drawCombo;
	GtkScale *specularity_scale;
} mlbrDialogWidgets;

static void query(void);
static void run(const gchar *name,
				gint nparams,
				const GimpParam *param,
				gint *nreturn_vals,
				GimpParam **return_vals);

static void mlBevelReflect(GimpDrawable *drawable,
				GimpPreview *preview);

static gboolean mlbr_dialog (GimpDrawable *drawable);

// Set up default values for options
static MLBevReflVals bVals = {
	0, // invert input
	3, // slope type
	3, // max normal map blur radius
	-1, // drawable id: invalid=normal map
	1, // specularity
	1  // preview
};

GimpPlugInInfo PLUG_IN_INFO = {
	NULL,
	NULL,
	query,
	run
};

static gboolean mlbr_constraint(gint32 imageId, gint32 drawableId, gpointer data);

MAIN()

static void
query (void) {
	static GimpParamDef args[] = {
		{GIMP_PDB_INT32,    "run-mode",         "The run mode { RUN-INTERACTIVE (0), RUN-NONINTERACTIVE (1) }"},
		{GIMP_PDB_IMAGE,    "image",            "Input image (unused)"},
		{GIMP_PDB_DRAWABLE, "drawable",         "Input and result drawable"},
		{GIMP_PDB_INT8,     "invertInput",      "Invert input"},
		{GIMP_PDB_INT8,     "slopeType",        "Slope type { flat (0), round (1)}"},
		{GIMP_PDB_INT8,     "maxNormMBlurRad",  "Maximum blur radius of depth map for flat slope"},
		{GIMP_PDB_DRAWABLE, "probeImgDrawable", "Spherical probe image drawable {none (-1)}"},
		{GIMP_PDB_FLOAT,    "specularity",      "Specularity of probe image reflection {diffuse (0) .. specular (1)}"}
	};

	gimp_install_procedure (
		PLUGIN_PROC,
		"ML Bevel Reflect",
		"Bevel white foreground object and apply reflection of probe image",
		"Marcel Lancelle",
		"Copyright Marcel Lancelle",
		PLUGIN_VERSION,
		"_ML Bevel Reflect...",
		"RGB*", //"RGB*, GRAY*",
		GIMP_PLUGIN,
		G_N_ELEMENTS (args), 0,
		args, NULL);

	gimp_plugin_menu_register (PLUGIN_PROC, "<Image>/Filters/Map");
}

static void
run (const gchar      *name,
	 gint              nparams,
	 const GimpParam  *param,
	 gint             *nreturn_vals,
	 GimpParam       **return_vals)
{
	static GimpParam  values[1];
	GimpPDBStatusType status = GIMP_PDB_SUCCESS;
	GimpRunMode       run_mode;
	GimpDrawable     *drawable;

	// Setting mandatory output values
	*nreturn_vals = 1;
	*return_vals  = values;

	values[0].type = GIMP_PDB_STATUS;
	values[0].data.d_status = status;

	// Getting run_mode - we won't display a dialog if 
	// we are in NONINTERACTIVE mode
	run_mode = param[0].data.d_int32;

	//  Get the specified drawable
	drawable = gimp_drawable_get(param[2].data.d_drawable);

	switch (run_mode) {
		case GIMP_RUN_INTERACTIVE:
			// Get options last values if needed
			gimp_get_data(PLUGIN_PROC, &bVals);

			// Display the dialog
			if (!mlbr_dialog(drawable)) return;
			break;

		case GIMP_RUN_NONINTERACTIVE:
			if (nparams != 8) {
				status = GIMP_PDB_CALLING_ERROR;
			} else {
				bVals.invertInput = param[3].data.d_int8;
				bVals.slopeType = param[4].data.d_int32;
				bVals.maxNormalMapBlurRadius = param[5].data.d_int32;
				bVals.probeImgDrawable = param[6].data.d_drawable;
				bVals.specularity = param[7].data.d_float;
			}
			break;

		case GIMP_RUN_WITH_LAST_VALS:
			gimp_get_data(PLUGIN_PROC, &bVals);
			break;

		default:
			break;
	}

	if (status == GIMP_PDB_SUCCESS) {
		mlBevelReflect (drawable, NULL);

		gimp_displays_flush ();
		gimp_drawable_detach (drawable);

		// Finally, set options in the core
		if (run_mode == GIMP_RUN_INTERACTIVE) {
			gimp_set_data(PLUGIN_PROC, &bVals, sizeof (MLBevReflVals));
		}
	}
}

static void
mlBevelReflect (GimpDrawable *drawable,
	  GimpPreview  *preview)
{
	// general temp variables
	gint i, x, y;
	guchar val;

	// temp vars for distance transformation
	gfloat h, v;
	gfloat hn1, vn1, sn1, hn2, vn2, sn2, hn3, vn3, sn3;
	gfloat vleft, vright, vtop, vbottom;

	gint channels, bytesPerRow;
	gboolean hasAlpha;
	gint x1, y1, x2, y2; // bounding rectangle for calculations
	GimpPixelRgn rgn_in, rgn_out;
	guchar *row, *rowAbove, *rowBelow, *tmpRow; // cache rows for input/output
	gint width, height;
	gfloat *distHoriz, *distVert, *distDiag; // for distance transform

	// temp vars for normal map
	gfloat dx, dy;
	gfloat numx, numy;
	gint kr, krsq, ksq, kx, ky;
	gfloat klen, df, ang;
	gfloat lx, ly;
	gint lix, liy;
	GimpDrawable *probeImgDrawable = NULL;
	GimpPixelRgn probeImgRgn;
	gint probeImgChannels;
	gboolean probeImgHasAlpha;
	guchar pix[4];

	gfloat nx, ny, nz, lnx, lny, lnz, tmpLen;
	gint inc;

	// specularity
	gint j;
	gint pixAvg[4];

	if (!preview) gimp_progress_init ("ML Bevel Reflect...");

	// get upper left and lower right coordinates
	if (preview) {
		gimp_preview_get_position(preview, &x1, &y1);
		gimp_preview_get_size(preview, &width, &height);
		x2 = x1 + width;
		y2 = y1 + height;
	} else {
		gimp_drawable_mask_bounds(drawable->drawable_id, &x1, &y1, &x2, &y2);
		width = x2 - x1;
		height = y2 - y1;
	}

	if ((!width) || (!height)) return; // can this happen anyway?

	if (!preview) gimp_progress_update(0.01);

	channels = gimp_drawable_bpp(drawable->drawable_id);
	hasAlpha = gimp_drawable_has_alpha(drawable->drawable_id);

	// allocate a big enough tile cache to read or write one row
	gimp_tile_cache_ntiles ((drawable->width / gimp_tile_width () + 1) + 16); // plus 16 for probe image , TODO

	// init two PixelRgns, one to read original data,
	// and the other to write output data
	gimp_pixel_rgn_init (&rgn_in, drawable,
						x1, y1, width, height,
						FALSE, FALSE);

	// allocate memory for input rows (current and the one above)
	bytesPerRow = width*channels;
	row = g_new(guchar, bytesPerRow);
	rowAbove = g_new(guchar, bytesPerRow);
	rowBelow = g_new(guchar, bytesPerRow);
	// allocate memory for distance transform
	distHoriz = g_new(gfloat, width*height);
	distVert = g_new(gfloat, width*height);

	if (!preview) gimp_progress_update(0.02);

	gimp_pixel_rgn_get_row(&rgn_in, row, x1, y1, width);
	gimp_pixel_rgn_get_row(&rgn_in, rowBelow, x1, y1+1, width);

	// initialize distance transform buffers
	for (y=0; y<height; y++) {
		// read current row
		tmpRow = rowAbove;
		rowAbove = row;
		row = rowBelow;
		rowBelow = tmpRow;
		gimp_pixel_rgn_get_row(&rgn_in, rowBelow, x1, y1 + MIN(y+1, height-1), width);

		for (x=0; x<width; x++) {
			// for now, use only first/red channel

			guchar val = row[x*channels];
			if (bVals.invertInput) val = 255-val;

			i = x + y*width;
			if (val == 255) {
				distHoriz[i] = 9999;
				distVert[i] = 9999;
			} else {
				distHoriz[i] = 0;
				distVert[i] = 0;
			}

			// trying to be smart with anti aliased edges
			if ((val > 0) && (val < 255)) {
				if ((x > 0) && (x < width-1) && (y > 0) && (y < height-1)) {
					vleft = row[(x-1)*channels];
					vright = row[(x+1)*channels];
					vtop = rowAbove[x*channels];
					vbottom = rowBelow[x*channels];
					if (MIN(vleft, vright) < MIN(vtop, vbottom)) {
						distHoriz[i] = val/255.0;
						distVert[i] = 0;
					} else {
						distHoriz[i] = 0;
						distVert[i] = val/255.0;
					}
				}
			}
		}
	}

	g_free(rowAbove);
	g_free(rowBelow);

	if (!preview) gimp_progress_set_text("ML Bevel Reflect: Distance Transform...");
	if (!preview) gimp_progress_update(0.05);

	// distance transform, first pass: top left to bottom right
	for (y=1; y<height; y++) {
		for (x=1; x<width-1; x++) {
			i = x + y*width;
			h = distHoriz[i];
			v = distVert[i];
			if ((h>0) || (v>0)) {
				i = x-1 + y*width;
				hn1 = distHoriz[i];
				hn1++;
				vn1 = distVert[i];
				sn1 = hn1*hn1 + vn1*vn1;

				i = x + (y-1)*width;
				hn2 = distHoriz[i];
				vn2 = distVert[i];
				vn2++;
				sn2 = hn2*hn2 + vn2*vn2;

				i = x+1 + (y-1)*width;
				hn3 = distHoriz[i];
				vn3 = distVert[i];
				hn3++;
				vn3++;
				sn3 = hn3*hn3 + vn3*vn3;

				if (sn2 < sn1) {
					sn1 = sn2;
					hn1 = hn2;
					vn1 = vn2;
				}

				if (sn3 < sn1) {
					sn1 = sn3;
					hn1 = hn3;
					vn1 = vn3;
				}

				if (sn1 < v*v + h*h) {
					i = x + y*width;
					distHoriz[i] = hn1;
					distVert[i] = vn1;
				}
			}
		}
		if ((!preview) && (y % 10 == 0)) gimp_progress_update((gdouble) y / (height-1) * 0.15 + 0.05);
	}

	// distance transform, second pass: bottom right to top left
	for (y=height-2; y>=0; y--) {
		for (x=width-2; x>0; x--) {
			i = x + y*width;
			h = distHoriz[i];
			v = distVert[i];
			if ((h>0) || (v>0)) {
				i = x+1 + y*width;
				hn1 = distHoriz[i];
				hn1++;
				vn1 = distVert[i];
				sn1 = hn1*hn1 + vn1*vn1;

				i = x + (y+1)*width;
				hn2 = distHoriz[i];
				vn2 = distVert[i];
				vn2++;
				sn2 = hn2*hn2 + vn2*vn2;

				i = x-1 + (y+1)*width;
				hn3 = distHoriz[i];
				vn3 = distVert[i];
				hn3++;
				vn3++;
				sn3 = hn3*hn3 + vn3*vn3;

				if (sn2 < sn1) {
					sn1 = sn2;
					hn1 = hn2;
					vn1 = vn2;
				}

				if (sn3 < sn1) {
					sn1 = sn3;
					hn1 = hn3;
					vn1 = vn3;
				}

				if (sn1 < v*v + h*h) {
					i = x + y*width;
					distHoriz[i] = hn1;
					distVert[i] = vn1;
				}
			}
		}
		if ((!preview) && (y % 10 == 0)) gimp_progress_update((gdouble) y / (height-1) * 0.15 + 0.20);
	}

	// distance transform, final pass
	distDiag = g_new(gfloat, width*height);
	for (y=0; y<height; y++) {
		for (x=0; x<width; x++) {
			i = x + y*width;
			h = distHoriz[i];
			v = distVert[i];
			distDiag[i] = sqrt(h*h + v*v);
		}
		if ((!preview) && (y % 10 == 0)) gimp_progress_update((gdouble) y / (height-1) * 0.15 + 0.35);
	}
	g_free(distHoriz);
	g_free(distVert);

	/////////////////

	/*
	// show distance transform
	gimp_pixel_rgn_init (&rgn_out,
					   drawable,
					   x1, y1,
					   width, height,
					   preview == NULL, TRUE);
	for (y=0; y<height; y++) {
		for (x=0; x<width; x++) {
			i = x + y*width;
			val = CLAMP((int) (20.f * distDiag[i]), 0, 255);
			row[x*channels] = val;
		}
		gimp_pixel_rgn_set_row (&rgn_out, row, x1, y1 + y, width);
	}
	*/

	/////////////////

	// distance transformation is now stored in distDiag
	// now interpret distDiag as height field

	if (!preview) gimp_progress_set_text("ML Bevel Reflect: Computing height field...");
	if (!preview) gimp_progress_update(0.50);

	// manipulate height field
	if (bVals.slopeType == SLOPE_ROUND) {
		for (y=0; y<height; y++)
			for (x=0; x<width; x++) {
				i = x + y*width;
				distDiag[i] = 7.f*sqrt(distDiag[i]);
			}
	}

	/////////////////

	if (!preview) gimp_progress_set_text("ML Bevel Reflect: Computing normal map...");
	if (!preview) gimp_progress_update(0.55);

	gimp_pixel_rgn_init(&rgn_out, drawable, x1, y1, width, height, preview == NULL, TRUE);

	//if (bVals.probeImgDrawable >= 0) {
	probeImgDrawable = NULL;
	if (gimp_drawable_is_valid(bVals.probeImgDrawable)) {
		probeImgDrawable = gimp_drawable_get(bVals.probeImgDrawable);
		gimp_pixel_rgn_init(&probeImgRgn, probeImgDrawable, 0, 0, probeImgDrawable->width, probeImgDrawable->height, FALSE, FALSE);
		probeImgChannels = gimp_drawable_bpp(probeImgDrawable->drawable_id);
		probeImgHasAlpha = gimp_drawable_has_alpha(probeImgDrawable->drawable_id);
	}

	// compute smoothed normal map
	for (y=1; y<height-1; y++) {
		for (x=1; x<width-1; x++) {
			i = x + y*width;
			if (distDiag[i] > 0) { // foreground
				if (bVals.slopeType == SLOPE_ROUND2) {
					// compute average gradient
					nx = 0;
					ny = 0;
					nz = 0;
					// radius for smoothing
					kr = (int) (distDiag[i])+2; // radius for circular region where to compute average gradient
					krsq = kr*kr;
					inc = (int) MAX(1, kr/6);
					kr -= kr%inc; // remove offset and thus artefacts

					for (ky=-kr; ky<=kr; ky++)
						for (kx=-kr; kx<=kr; kx++) { // for all pixels in that square region
							if ((x+kx >= 0) && (y+ky >= 0) && (x+kx < width) && (y+ky < height)) {
								ksq = kx*kx + ky*ky;
								if ((ksq <= krsq) && ((kx != 0) || (ky != 0))) { // inside circle (larger?) and not 0,0
									klen = sqrt((float) ksq);
									df = 20.f*(sqrt(distDiag[x+kx + (y+ky)*width]) - sqrt(distDiag[i]));
									lnx = -df*kx/klen;
									lny = -df*ky/klen;
									lnz = klen;
									tmpLen = sqrt(lnx*lnx + lny*lny + lnz*lnz);
									nx += lnx/tmpLen;
									ny += lny/tmpLen;
									nz += lnz/tmpLen;
								}
							}
						}
					// normalize
					tmpLen = sqrt(nx*nx + ny*ny + nz*nz);
					lx = nx/tmpLen;
					ly = ny/tmpLen;
				} else {
					// compute average gradient
					dx = 0;
					dy = 0;
					numx = 0.01;
					numy = 0.01;
					// radius for smoothing
					kr = (int) (distDiag[i])+2; // radius for circular region where to compute average gradient
					if (bVals.slopeType == SLOPE_FLAT) {
						if (kr > bVals.maxNormalMapBlurRadius) kr = bVals.maxNormalMapBlurRadius;
					}
					krsq = kr*kr;
					for (ky=-kr; ky<=kr; ky++)
						for (kx=-kr; kx<=kr; kx++) { // for all pixels in that square region
							if ((x+kx >= 0) && (y+ky >= 0) && (x+kx < width) && (y+ky < height)) {
								ksq = kx*kx + ky*ky;
								if ((ksq <= krsq) && ((kx != 0) || (ky != 0))) { // inside circle (larger?) and not 0,0
									klen = sqrt((float) ksq);
									df = (distDiag[x+kx + (y+ky)*width] - distDiag[i]) / klen;
									//if (distDiag[y+ky][x+kx] == 0) df *= 4; //df = (df > 0) ? 4 : -4;
									ang = atan2((float) ky, (float) kx);
									dx += cos(ang) * df;
									dy += sin(ang) * df;
									numx += fabs(cos(ang));
									numy += fabs(sin(ang));
								}
							}
						}

					lx = -cos(atan2(1,dx/numx));
					ly = -cos(atan2(1,dy/numy));
				}

				i = x*channels;
				if (probeImgDrawable) {
						// experimental ////////
						if (bVals.specularity < 1) {
							kr = (int) MIN(probeImgDrawable->width*(1-bVals.specularity*bVals.specularity)*cos(0.5*3.1415926*lx),
											probeImgDrawable->height*(1-bVals.specularity*bVals.specularity)*cos(0.5*3.1415926*ly));
											lix = (int) (0.5*(probeImgDrawable->width-1)*(1+0.99*lx));
											liy = (int) (0.5*(probeImgDrawable->height-1)*(1+0.99*ly));
						pixAvg[0]=0;
						pixAvg[1]=0;
						pixAvg[2]=0;
						pixAvg[3]=0;
						for (j=0; j<8; j++) {
							gimp_pixel_rgn_get_pixel(&probeImgRgn, pix, lix + kr*(rand()/(float)RAND_MAX-0.5), liy + kr*(rand()/(float)RAND_MAX-0.5));
							pixAvg[0]+=pix[0];
							pixAvg[1]+=pix[1];
							pixAvg[2]+=pix[2];
							pixAvg[3]+=pix[3];
						}
						row[i++] = pixAvg[0]/8;
						row[i++] = pixAvg[1]/8;
						row[i++] = pixAvg[2]/8;
						row[i] = probeImgHasAlpha ? pixAvg[3]/8 : 255;
						/////////////////////////
					} else {
						lix = (int) (0.5*(probeImgDrawable->width-1)*(1+0.99*lx));
						liy = (int) (0.5*(probeImgDrawable->height-1)*(1+0.99*ly));
						gimp_pixel_rgn_get_pixel(&probeImgRgn, pix, lix, liy);
						row[i++] = pix[0];
						row[i++] = pix[1];
						row[i++] = pix[2];
						if (hasAlpha) row[i] = probeImgHasAlpha ? pix[3] : 255;
					}
				} else {
					// normal map
					row[i++] = 127 + (int) (127*lx);
					row[i++] = 127 + (int) (127*ly);
					row[i++] = 127 + (int) (127*(1-sqrt(lx*lx+ly*ly)));
					if (hasAlpha) row[i] = 255;
				}
			} else { // background
				i = x*channels;
				if (probeImgDrawable) {
					// use color of top left corner of probe image for background
					gimp_pixel_rgn_get_pixel(&probeImgRgn, pix, 0, 0);
					row[i++] = pix[0];
					row[i++] = pix[1];
					row[i++] = pix[2];
					if (hasAlpha) row[i] = probeImgHasAlpha ? pix[3] : 255;
				} else {
					row[i++] = 127;
					row[i++] = 127;
					row[i++] = 255;
					if (hasAlpha) row[i] = 255;
				}
			}
		}
		gimp_pixel_rgn_set_row (&rgn_out, row, x1, y1 + y, width);

		if ((!preview) && (y % 10 == 0)) gimp_progress_update((gdouble) y / (height-1) * 0.4 + 0.59);
	}

	// gimp_progress_end(); not needed as plugin terminates anyway

	if (probeImgDrawable != NULL) gimp_drawable_detach(probeImgDrawable);

	// free memory
	g_free(distDiag);
	g_free(row);

	//  Update the modified region
	if (preview) {
		gimp_drawable_preview_draw_region(GIMP_DRAWABLE_PREVIEW (preview), &rgn_out);
	} else {
		gimp_drawable_flush(drawable);
		gimp_drawable_merge_shadow(drawable->drawable_id, TRUE);
		gimp_drawable_update(drawable->drawable_id, x1, y1, width, height);
	}
}

// callbacks ///////////////////////////////////////////

static gboolean
mlpr_constraint (gint32 imageId,
				gint32 drawableId,
				gpointer data)
{
	return ((drawableId == -1) || gimp_drawable_is_rgb(drawableId) );
}

void invertInput_changed_callback(GtkCheckButton *checkbutton, gpointer data) {
	bVals.invertInput = gtk_toggle_button_get_active(mlbrDialogWidgets.invertInput_checkButton);
	gimp_preview_invalidate(mlbrDialogWidgets.preview);
}

void slopeType_changed_callback(GtkRadioButton *radiobutton, gpointer data) { 
	if (gtk_toggle_button_get_active(mlbrDialogWidgets.slope_flat_radiobutton)) {
		bVals.slopeType = SLOPE_FLAT;
		gimp_preview_invalidate(mlbrDialogWidgets.preview);
	}
	if (gtk_toggle_button_get_active(mlbrDialogWidgets.slope_round_radiobutton)) {
		bVals.slopeType = SLOPE_ROUND;
		gimp_preview_invalidate(mlbrDialogWidgets.preview);
	}
	if (gtk_toggle_button_get_active(mlbrDialogWidgets.slope_round2_radiobutton)) {
		bVals.slopeType = SLOPE_ROUND2;
		gimp_preview_invalidate(mlbrDialogWidgets.preview);
	}
}

void reflection_changed_callback(GtkRadioButton *unused, gpointer data) {
	if (gtk_toggle_button_get_active(mlbrDialogWidgets.probeImg_radiobutton)) {
		//bVals.probeImgDrawable = gtk_combo_box_get_active(mlbrDialogWidgets.probeImg_drawCombo);
		gimp_int_combo_box_get_active(GIMP_INT_COMBO_BOX(mlbrDialogWidgets.probeImg_drawCombo), &bVals.probeImgDrawable);
		gimp_preview_invalidate(mlbrDialogWidgets.preview);
	}
	if (gtk_toggle_button_get_active(mlbrDialogWidgets.normalMap_radiobutton)) {
		bVals.probeImgDrawable = -1;
		gimp_preview_invalidate(mlbrDialogWidgets.preview);
	}
}

void probeImgDrawableChangedCallback(GtkWidget *widget, gpointer preview) {
	//g_message("probeImgDrawab Callback\n");
	//gimp_preview_invalidate(mlbrDialogWidgets.preview);
	reflection_changed_callback(NULL, NULL);
}

void specularity_value_changed(GtkWidget* scale, gpointer data) {
	bVals.specularity = gtk_range_get_value(GTK_RANGE(mlbrDialogWidgets.specularity_scale));
	gimp_preview_invalidate(mlbrDialogWidgets.preview);
}

// dialog /////////////////////////////////////////////

static gboolean
mlbr_dialog (GimpDrawable *drawable) {
	GtkWidget *dialog;
	GtkWidget *main_vbox;
		GtkWidget *preview;
		GtkWidget *main_hbox;
			GtkWidget *left_vbox;
				GtkWidget *input_frame;
					GtkWidget *invertInput_checkButton;
				GtkWidget *slope_frame;
					GtkWidget *slope_vbox;
						GtkWidget *slope_flat_hbox;
							GtkWidget *slope_flat_radiobutton;
								GtkWidget *flat_frame;
									GtkWidget *flat_hbox;
										GtkWidget *flat_radius_spinbutton;
										GtkObject *flat_radius_spinbutton_adj;
						GtkWidget *slope_round_radiobutton;
						GtkWidget *slope_round2_radiobutton;
			GtkWidget *reflection_frame;
				GtkWidget *reflection_vbox;
					GtkWidget *normalMap_radiobutton;
					GtkWidget *ref_probeImg_hbox;
						GtkWidget *probeImg_radiobutton;
							GtkWidget *probeImg_vbox;
								GtkWidget *probeImg_drawCombo;
								GtkWidget *specularity_frame;
									GtkWidget *specularity_vbox;
										GtkWidget *specularity_scale;
										GtkObject *specularity_scale_adj;
										GtkWidget *specularity_hbox;

	GtkWidget *tmp_label;
	//GtkWidget *tmp_separator;
	GtkWidget * tmp_alignment;

	gboolean run;

	gimp_ui_init ("mlbevelreflect", FALSE);

	dialog = gimp_dialog_new ("ML Bevel Reflect", "mlbevelreflect",
							NULL, 0,
							gimp_standard_help_func, PLUGIN_PROC, // TODO: provide help html site??

							GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
							GTK_STOCK_OK,     GTK_RESPONSE_OK,

							NULL);

	main_vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 6);
	gtk_container_add(GTK_CONTAINER (GTK_DIALOG(dialog)->vbox), main_vbox);

	preview = gimp_drawable_preview_new(drawable, &bVals.preview);
	mlbrDialogWidgets.preview = (GimpPreview*) preview;
	gtk_box_pack_start(GTK_BOX(main_vbox), preview, TRUE, TRUE, 0);
	g_signal_connect_swapped (preview, "invalidated", G_CALLBACK (mlBevelReflect), drawable);

	main_hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
	gtk_box_pack_start(GTK_BOX(main_vbox), main_hbox, FALSE, FALSE, 0);

	left_vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 2);
	gtk_box_pack_start(GTK_BOX(main_hbox), left_vbox, FALSE, FALSE, 0);

	input_frame = gtk_frame_new(NULL);
	gtk_box_pack_start(GTK_BOX(left_vbox), input_frame, TRUE, TRUE, 2);
	gtk_container_set_border_width(GTK_CONTAINER(input_frame), 6);
	tmp_label = gtk_label_new(NULL);
	gtk_label_set_markup((GtkLabel*) tmp_label, "<b>Input</b>");
	gtk_frame_set_label_widget(GTK_FRAME(input_frame), tmp_label);

	invertInput_checkButton = gtk_check_button_new_with_mnemonic("_Invert input");
	mlbrDialogWidgets.invertInput_checkButton = (GtkCheckButton*) invertInput_checkButton;
	gtk_container_add(GTK_CONTAINER(input_frame), invertInput_checkButton);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(invertInput_checkButton), bVals.invertInput);
	g_signal_connect(invertInput_checkButton, "toggled", G_CALLBACK(invertInput_changed_callback), NULL);

	slope_frame = gtk_frame_new(NULL);
	gtk_box_pack_start(GTK_BOX(left_vbox), slope_frame, TRUE, TRUE, 2);
	gtk_container_set_border_width(GTK_CONTAINER(slope_frame), 6);
	tmp_label = gtk_label_new(NULL);
	gtk_label_set_markup((GtkLabel*) tmp_label, "<b>Slope Options</b>");
	gtk_frame_set_label_widget(GTK_FRAME(slope_frame), tmp_label);

	slope_vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 2);
	gtk_container_add(GTK_CONTAINER(slope_frame), slope_vbox);

	slope_flat_hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	gtk_box_pack_start(GTK_BOX(slope_vbox), slope_flat_hbox, TRUE, TRUE, 0);

	slope_flat_radiobutton = gtk_radio_button_new_with_mnemonic(NULL, "_Flat");
	mlbrDialogWidgets.slope_flat_radiobutton = (GtkRadioButton*) slope_flat_radiobutton;
	gtk_box_pack_start(GTK_BOX(slope_flat_hbox), slope_flat_radiobutton, TRUE, TRUE, 0);

	if (bVals.slopeType == SLOPE_FLAT) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(slope_flat_radiobutton), TRUE);
	g_signal_connect(slope_flat_radiobutton, "toggled", G_CALLBACK(slopeType_changed_callback), NULL);

	flat_frame = gtk_frame_new(NULL);
	gtk_box_pack_start(GTK_BOX(slope_flat_hbox), flat_frame, TRUE, TRUE, 2);

	flat_hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
	gtk_container_add(GTK_CONTAINER(flat_frame), flat_hbox);

	tmp_label = gtk_label_new_with_mnemonic("_Max. smoothing radius");
	gtk_box_pack_start(GTK_BOX(flat_hbox), tmp_label, FALSE, FALSE, 0);
	flat_radius_spinbutton = gimp_spin_button_new(&flat_radius_spinbutton_adj,
		bVals.maxNormalMapBlurRadius, 1, 32, 1, 1, 1, 5, 0);
	gtk_box_pack_start(GTK_BOX(flat_hbox), flat_radius_spinbutton, FALSE, FALSE, 0);
	g_signal_connect_swapped(flat_radius_spinbutton_adj, "value_changed",
							G_CALLBACK (gimp_preview_invalidate), preview);
	g_signal_connect(flat_radius_spinbutton_adj, "value_changed",
					G_CALLBACK(gimp_int_adjustment_update), &bVals.maxNormalMapBlurRadius);

	slope_round_radiobutton = gtk_radio_button_new_with_mnemonic_from_widget(GTK_RADIO_BUTTON(slope_flat_radiobutton), "_Round");
	mlbrDialogWidgets.slope_round_radiobutton = (GtkRadioButton*) slope_round_radiobutton;
	gtk_box_pack_start(GTK_BOX(slope_vbox), slope_round_radiobutton, TRUE, TRUE, 0);
	if (bVals.slopeType == SLOPE_ROUND) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(slope_round_radiobutton), TRUE);
	g_signal_connect(slope_round_radiobutton, "toggled", G_CALLBACK(slopeType_changed_callback), NULL);

	slope_round2_radiobutton = gtk_radio_button_new_with_mnemonic_from_widget(GTK_RADIO_BUTTON(slope_flat_radiobutton), "_Round v2");
	mlbrDialogWidgets.slope_round2_radiobutton = (GtkRadioButton*) slope_round2_radiobutton;
	gtk_box_pack_start(GTK_BOX(slope_vbox), slope_round2_radiobutton, TRUE, TRUE, 0);
	if (bVals.slopeType == SLOPE_ROUND2) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(slope_round2_radiobutton), TRUE);
	g_signal_connect(slope_round2_radiobutton, "toggled", G_CALLBACK(slopeType_changed_callback), NULL);

	reflection_frame = gtk_frame_new(NULL);
	gtk_box_pack_start(GTK_BOX(main_hbox), reflection_frame, TRUE, TRUE, 2);
	gtk_container_set_border_width(GTK_CONTAINER(reflection_frame), 6);
	tmp_label = gtk_label_new(NULL);
	gtk_label_set_markup((GtkLabel*) tmp_label, "<b>Reflection Color</b>");
	gtk_frame_set_label_widget(GTK_FRAME(reflection_frame), tmp_label);

	reflection_vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 2);
	gtk_container_add(GTK_CONTAINER(reflection_frame), reflection_vbox);

	normalMap_radiobutton = gtk_radio_button_new_with_mnemonic(NULL, "_Normal map");
	mlbrDialogWidgets.normalMap_radiobutton = (GtkRadioButton*) normalMap_radiobutton;
	gtk_box_pack_start(GTK_BOX(reflection_vbox), normalMap_radiobutton, TRUE, TRUE, 0);
	if (bVals.probeImgDrawable < 0) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(normalMap_radiobutton), TRUE);
	g_signal_connect(normalMap_radiobutton, "toggled", G_CALLBACK(reflection_changed_callback), NULL);

	ref_probeImg_hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 2);
	gtk_box_pack_start(GTK_BOX(reflection_vbox), ref_probeImg_hbox, TRUE, TRUE, 0);

	tmp_alignment = gtk_alignment_new(0, 0, 0, 0); // top left
	gtk_box_pack_start(GTK_BOX(ref_probeImg_hbox), tmp_alignment, FALSE, FALSE, 0);

	probeImg_radiobutton = gtk_radio_button_new_with_mnemonic_from_widget(GTK_RADIO_BUTTON(normalMap_radiobutton), "_Spherical probe image");
	mlbrDialogWidgets.probeImg_radiobutton = (GtkRadioButton*) probeImg_radiobutton;
	gtk_container_add(GTK_CONTAINER(tmp_alignment), probeImg_radiobutton);
	if (bVals.probeImgDrawable >= 0) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(probeImg_radiobutton), TRUE);
	g_signal_connect(probeImg_radiobutton, "toggled", G_CALLBACK(reflection_changed_callback), NULL);

	probeImg_vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
	gtk_box_pack_start(GTK_BOX(ref_probeImg_hbox), probeImg_vbox, TRUE, TRUE, 0);

	probeImg_drawCombo = gimp_drawable_combo_box_new(mlpr_constraint, NULL);
	mlbrDialogWidgets.probeImg_drawCombo = (GimpDrawableComboBox*) probeImg_drawCombo;
	gtk_box_pack_start(GTK_BOX(probeImg_vbox), probeImg_drawCombo, FALSE, FALSE, 0);
	gimp_int_combo_box_connect(GIMP_INT_COMBO_BOX(probeImg_drawCombo), bVals.probeImgDrawable,
								G_CALLBACK(probeImgDrawableChangedCallback), NULL);

	specularity_frame = gtk_frame_new(NULL);
	gtk_box_pack_start(GTK_BOX(probeImg_vbox), specularity_frame, TRUE, TRUE, 0);
	gtk_container_set_border_width(GTK_CONTAINER(specularity_frame), 6);
	tmp_label = gtk_label_new(NULL);
	gtk_label_set_markup((GtkLabel*) tmp_label, "<b>Specularity</b>");
	gtk_frame_set_label_widget(GTK_FRAME(specularity_frame), tmp_label);

	specularity_vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 2);
	gtk_container_add(GTK_CONTAINER(specularity_frame), specularity_vbox);

	specularity_scale_adj = gtk_adjustment_new(0.0, 0.0, 1.0, 0.1, 0.5, 0.0);
	specularity_scale = gtk_hscale_new(GTK_ADJUSTMENT(specularity_scale_adj));
	mlbrDialogWidgets.specularity_scale = (GtkScale*) specularity_scale;
	gtk_range_set_value(GTK_RANGE(mlbrDialogWidgets.specularity_scale), bVals.specularity);
	gtk_box_pack_start(GTK_BOX(specularity_vbox), specularity_scale, TRUE, TRUE, 0);
	g_signal_connect(specularity_scale, "value_changed",
		G_CALLBACK(specularity_value_changed), NULL);

	specularity_hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 2);
	gtk_box_pack_start(GTK_BOX(specularity_vbox), specularity_hbox, TRUE, TRUE, 0);

	tmp_label = gtk_label_new("Diffuse");
	gtk_box_pack_start(GTK_BOX(specularity_hbox), tmp_label, FALSE, FALSE, 0);

	tmp_label = gtk_label_new("Shiny");
	gtk_box_pack_end(GTK_BOX(specularity_hbox), tmp_label, FALSE, FALSE, 0);


	mlBevelReflect(drawable, GIMP_PREVIEW(preview));

	gtk_widget_show_all (dialog);

	run = (gimp_dialog_run(GIMP_DIALOG(dialog)) == GTK_RESPONSE_OK);

	gtk_widget_destroy(dialog);

	return run;
}
