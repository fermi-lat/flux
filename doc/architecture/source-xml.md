# flux source XML reference

This page describes the source-library XML that `FluxMgr` reads, as the
code in `src/FluxMgr.cxx` and `src/FluxSource.cxx` interprets it at
`1d5fc58`. `xml/source.dtd` is the intended grammar, but it is **not read**
any more (RapidXML does no DTD processing), so its default values do not
apply. See [README.md](README.md) for the overall design.

## Structure

```xml
<source_library title="my sources">          <!-- title is reported by sourceOriginList() -->

  <!-- 1. A simple particle source: flux on <source>, spectrum built in -->
  <source name="pt_src" flux="1e-3" ident="7">
    <spectrum escale="MeV" apply_edisp="true">
      <particle name="gamma">
        <power_law emin="30" emax="2e5" gamma="2.1"/>
      </particle>
      <celestial_dir ra="83.6" dec="22.0"/>
    </spectrum>
  </source>

  <!-- 2. A plug-in spectrum class with a free-form params string -->
  <source name="crab_pulsar">
    <spectrum escale="MeV">
      <SpectrumClass name="Pulsar"
                     params="1e-3,2.,0.033,0,0,$(GENERICSOURCESDATAPATH)/CrabTemplate.dat"/>
      <celestial_dir ra="83.57" dec="22.01"/>
    </spectrum>
  </source>

  <!-- 3. A spectrum class that also supplies directions -->
  <source name="gal_diffuse">
    <spectrum escale="MeV">
      <SpectrumClass name="MapCube" params="flux=17.,fitsFile=$(MYDATA)/gll_iem.fits"/>
      <use_spectrum frame="galaxy"/>
    </spectrum>
  </source>

  <!-- 4. A composite (see the gotcha: only the first nestedSource is used today) -->
  <source name="both">
    <nestedSource sourceRef="pt_src"/>
    <nestedSource sourceRef="crab_pulsar"/>
  </source>
</source_library>
```

## Elements and attributes

The "DTD default" column is what the DTD declares; the "Code today" column
is what the code uses when the attribute is missing.

| Element | Attribute | Meaning | DTD default | Code today |
|---------|-----------|---------|-------------|------------|
| `source_library` | `title` | Label for `sourceOriginList()` | `user` | file name |
| `source` | `name` | Lookup key (`FluxMgr::source(name)`) | required | — |
| `source` | `flux` | Flux used by `<particle>` sources (particles/m²/s/sr) | `1.0` | `0` (and ignored, see gotcha) |
| `source` | `ident` | Positive value becomes the spectrum `identifier()` → `MC_SRC_ID` | `-1` | `0` |
| `spectrum` | `escale` | `MeV` or `GeV`. Anything else aborts construction | required | — |
| `spectrum` | `flux` | Passed to `Spectrum::setFlux` for `SpectrumClass` sources (only meaningful if the class uses `m_flux`, e.g. `FileSpectrum`) | none | none |
| `spectrum` | `apply_edisp` | `false` tells `gtobssim` not to apply energy dispersion to this source. `true`/empty → apply. Other values throw | `true` | apply |
| `spectrum` | `particle_name` | Particle name for `SpectrumClass` sources | `gamma` | **ignored** |
| `particle` | `name` | Particle name for `SimpleSpectrum` (`gamma`, `mu+`, …) | `gamma` | `""` |
| `energy` | `e` | Monoenergetic, in `escale` units | `0.1` | `0` |
| `power_law` | `emin`, `emax`, `gamma`, `ebreak`, `gamma2` | Power law; broken if `ebreak > 0` | `0.1, 100, 2, 0, 0` | all `0` |
| `SpectrumClass` | `name` | Registered class name | required | — |
| `SpectrumClass` | `params` | Passed verbatim to `T(const std::string&)`. `$(VAR)` is expanded **by the class**, not by flux | `" "` | `""` |
| `<AnyName/>` | | Shorthand for `<SpectrumClass name="AnyName" params=""/>` | | |

**Direction** (first match among the children of `<spectrum>`; if none, the
source is isotropic over the full sphere in the instrument frame):

| Element | Attributes | Behaviour | Occultable |
|---------|-----------|-----------|------------|
| `celestial_dir` | `ra`, `dec`, `radius` (deg) | Fixed sky direction; uniform disk if `radius>0` | yes |
| `galactic_dir` | `l`, `b`, `radius` (deg) | As above in Galactic coordinates | yes |
| `use_spectrum` | `frame` = `zenith` (default), `nadir`, `equatorial`, `galactic`/`galaxy`, `Sun`, `Moon`, `Jupiter`, `Saturn` | Direction from `ISpectrum::dir(energy)`: `(cosθ, φ)` for zenith/nadir, `(lon, lat)` deg for sky frames, offsets for solar-system bodies. Solid angle from `ISpectrum::solidAngle()` | yes unless `zenith` |
| `direction` | `theta`, `phi` (deg), `frame` (`zenith` or anything else = spacecraft) | Fixed direction | only if `zenith` |
| `solid_angle` | `mincos`, `maxcos`, `theta`, `phi` | Random in a band of zenith angle, optionally tilted | no |
| `custom_dir` | | Uses `ISpectrum::launchDirection()` | no |

**Launch point** (the element after the direction element; default is a
random point on a disk of area `EventSource::totalArea()`, backed off
`s_backoff` = 2000 mm along the direction):

| Element | Attributes (mm) |
|---------|-----------------|
| `launch_point` | `x`, `y`, `z`, `beam_radius` |
| `patch` | `xmin`, `xmax`, `ymin`, `ymax`, `zmin`, `zmax` |
| `custom_pt` | uses `ISpectrum::launchPoint()` |

`gtobssim` ignores launch points; they matter only for detector simulation.

**Composites:** a `<source>` whose first child element is
`<nestedSource sourceRef="X"/>` is resolved through the library by name.

## How a definition becomes photons

1. `FluxMgr::source(name)` → `getSourceFromXML(node)` →
   `new FluxSource(node)`.
2. `FluxSource` builds the spectrum: `SimpleSpectrum` for `<particle>`,
   otherwise `SpectrumFactoryTable::instance()->instantiate(name, params)`.
   An unknown name prints the registered names and throws.
3. It builds the direction strategy and the launch strategy.
4. Each `event(t)`:
   - `interval = spectrum->interval(t)`. If not positive,
     `interval = -ln(1-u) / (flux(t) · solidAngle · totalArea)`.
   - If `t + interval < GPS::endTime()`: `energy = spectrum->energy(t+Δt)`;
     the direction strategy computes the LAT-frame direction at `t+Δt`
     (via `astro::GPS`); the launch strategy computes the point.
   - Otherwise the source disables itself.
5. A `CompositeSource` keeps each child's next event in a time-ordered
   multimap and hands out the earliest one.

## Behaviour changes in the RapidXML rewrite

Commit `1d5fc58` (PR #4 "xerces-replacement", merged 2026-07-24) rewrote
`FluxMgr.cxx`, `FluxSource.cxx` and `SimpleSpectrum.cxx` for RapidXML. The
last Xerces version is `a069d66`. Comparing the two by reading the code
(nothing below was checked at runtime):

| Area | Before (`a069d66`) | Now (`1d5fc58`) | Effect on `gtobssim` |
|------|--------------------|-----------------|----------------------|
| Parsing | All files wrapped in one generated document with external entities, validated against `source.dtd`; parse error → `FATAL_MACRO` (throws) | Each file parsed separately, no DTD; errors print a warning and skip the file | Bad XML no longer stops the tool at load time |
| DTD defaults | Applied (e.g. `<source flux>`=1.0, `mincos`=−0.4) | Not applied; missing → `""`/`0` | Sources relying on defaults change |
| `nestedSource` | `CompositeSource` of every `nestedSource` child | Only the first `sourceRef` | Composite catalogue entries emit one member |
| `FluxSource::flux()` | Spectrum flux if non-zero, else the `<source flux>` attribute; 0 if disabled | Spectrum flux always | `<particle>` sources have rate 0 and never emit |
| `escale="GeV"` | Energy × 1000 in `computeLaunch` | No conversion | GeV-scale sources' energies off by 1000 |
| `name()` | `spectrum()->name()` if non-empty, else XML name | XML name | `SourcePopulation` members lose their names |
| Particle name | `<spectrum particle_name>` | `<particle name>` child only | `SpectrumClass` sources can't set a particle name |
| Filter cone | `occulted()` applied `s_cone` | Cone test removed | `setFilterCone` is a no-op |
| Misalignment | `s_alignMatrix` applied to the launch direction | Not applied | `setAlignmentRotation(…, true)` is a no-op |
| No direction element | Error ("Unknown angle specification") or a null dereference | Isotropic `RandomDirection(-1, 1)` | Silent default |
| `disable()` | Deleted the spectrum | Keeps it | Fewer dangling-pointer risks |
| `skyDirection()` | `SkyDir(lat_dir)` (wrong frame) | `GPS::toSky(lat_dir)` | Fix; not used by `gtobssim` |

Branch `origin/f26-cpp17-compliance` (2026-10-16) only changes function
signatures to trailing-return style; it does not address any of these.
