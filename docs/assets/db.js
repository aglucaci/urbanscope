(function(){
  const MANIFEST_URL = "./db/srr_records_manifest.json";

  function norm(value){
    return (value ?? "").toString().trim();
  }

  function safeInt(value){
    const n = parseInt((value ?? "").toString(), 10);
    return Number.isFinite(n) ? n : null;
  }

  function normalizePartPath(path){
    let s = (path ?? "").toString().replaceAll("\\", "/");
    if (s.startsWith("docs/")) s = s.slice(5);
    if (!s.includes("/")) s = "db/" + s;
    if (s.startsWith("db/")) s = "./" + s;
    if (!s.startsWith("./")) s = "./" + s.replace(/^\/+/, "");
    return s;
  }

  function yearFromRecord(rec){
    const raw = rec?.runinfo_row?.ReleaseDate || rec?.runinfo_row?.LoadDate || "";
    const match = raw.match(/^(\d{4})-/);
    return match ? safeInt(match[1]) : null;
  }

  function getRun(rec){
    return rec?.srr || rec?.runinfo_row?.Run || "";
  }

  function getBioProject(rec){
    return rec?.runinfo_row?.BioProject || rec?.bioproject?.accession || "";
  }

  function getBioSample(rec){
    return rec?.runinfo_row?.BioSample || rec?.geo?.biosample_accession || "";
  }

  function getCountry(rec){
    return rec?.geo?.country || "";
  }

  function getCity(rec){
    return rec?.geo?.city || "";
  }

  function getAssay(rec){
    return rec?.assay?.assay_class || rec?.runinfo_row?.LibraryStrategy || "Unknown";
  }

  function getCenter(rec){
    return rec?.runinfo_row?.CenterName || rec?.bioproject?.center_name || "";
  }

  function getTitle(rec){
    return rec?.bioproject?.title || rec?.title || "";
  }

  function isKnownGeo(value){
    const s = norm(value).toLowerCase();
    return !!s && s !== "(unknown)";
  }

  async function fetchManifest(){
    const res = await fetch(MANIFEST_URL, { cache: "no-store" });
    if (!res.ok) throw new Error("Unable to load manifest");
    return await res.json();
  }

  async function fetchAllRuns(progressCb){
    const manifest = await fetchManifest();
    const parts = manifest.parts || [];
    const allRuns = [];

    for (let i = 0; i < parts.length; i += 1){
      const url = normalizePartPath(parts[i].path);
      const res = await fetch(url, { cache: "no-store" });
      if (!res.ok) throw new Error("Unable to load " + url);
      const chunk = await res.json();
      if (Array.isArray(chunk)) allRuns.push(...chunk);
      if (progressCb) progressCb({ done: i + 1, total: parts.length, manifest, allRuns });
    }

    return { manifest, allRuns };
  }

  function groupProjects(allRuns){
    const map = new Map();
    for (const rec of allRuns){
      const bp = getBioProject(rec) || "(unassigned)";
      if (!map.has(bp)){
        map.set(bp, {
          accession: bp,
          title: getTitle(rec),
          records: [],
          runs: new Set(),
          biosamples: new Set(),
          countries: new Set(),
          cities: new Set(),
          assays: new Set(),
          centers: new Set(),
          years: new Set()
        });
      }
      const row = map.get(bp);
      row.records.push(rec);
      if (getRun(rec)) row.runs.add(getRun(rec));
      if (getBioSample(rec)) row.biosamples.add(getBioSample(rec));
      if (getCountry(rec)) row.countries.add(getCountry(rec));
      if (getCity(rec)) row.cities.add(getCity(rec));
      if (getAssay(rec)) row.assays.add(getAssay(rec));
      if (getCenter(rec)) row.centers.add(getCenter(rec));
      const year = yearFromRecord(rec);
      if (year) row.years.add(year);
      if (!row.title && getTitle(rec)) row.title = getTitle(rec);
    }
    return Array.from(map.values());
  }

  function tally(values){
    const counts = new Map();
    for (const value of values){
      const key = norm(value) || "(unknown)";
      counts.set(key, (counts.get(key) || 0) + 1);
    }
    return Array.from(counts.entries())
      .map(([name, count]) => ({ name, count }))
      .sort((a, b) => b.count - a.count || a.name.localeCompare(b.name));
  }

  function summarize(allRuns){
    const projects = groupProjects(allRuns);
    const years = tally(allRuns.map(yearFromRecord).filter(Boolean));
    const assays = tally(allRuns.map(getAssay));
    const countries = tally(allRuns.map(getCountry));
    const cities = tally(allRuns.map(getCity));
    const centers = tally(allRuns.map(getCenter));

    const geoResolvedRuns = allRuns.filter((rec) => isKnownGeo(getCountry(rec)) && isKnownGeo(getCity(rec))).length;
    const downloadableRuns = allRuns.filter((rec) => norm(rec?.runinfo_row?.download_path)).length;

    return {
      totalRuns: allRuns.length,
      totalProjects: projects.length,
      totalBioSamples: new Set(allRuns.map(getBioSample).filter(Boolean)).size,
      totalCountries: countries.filter((d) => isKnownGeo(d.name)).length,
      totalCities: cities.filter((d) => isKnownGeo(d.name)).length,
      geoResolvedRuns,
      downloadableRuns,
      years,
      assays,
      countries,
      cities,
      centers,
      projects
    };
  }

  function topItems(items, limit){
    return items.slice(0, limit);
  }

  function formatNumber(value){
    return new Intl.NumberFormat("en-US").format(value || 0);
  }

  window.UrbanScopeDB = {
    fetchManifest,
    fetchAllRuns,
    summarize,
    topItems,
    formatNumber,
    yearFromRecord,
    getCountry,
    getCity,
    getAssay,
    getCenter,
    getBioProject,
    getBioSample,
    groupProjects
  };
})();
