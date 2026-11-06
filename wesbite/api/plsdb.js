
const PLSDB_API =
  "https://ccb-microbe.cs.uni-saarland.de/plsdb2025/api/sequence/";

export default async function handler(req, res) {
  // Only accept submission and polling requests.
  if (req.method !== "POST" && req.method !== "GET") {
    res.setHeader("Allow", "GET, POST");
    return res.status(405).json({ error: "Method not allowed" });
  }

  try {
    const upstreamUrl = new URL(PLSDB_API);

    // Forward polling parameters, e.g. ?job_id=...
    if (req.method === "GET") {
      upstreamUrl.search = new URL(req.url, "http://localhost").search;
    }

    const headers = {
      Accept: "application/json"
    };

    // Preserve the multipart boundary when forwarding the FASTA upload.
    if (req.method === "POST") {
      const contentType = req.headers["content-type"];

      if (!contentType || !contentType.includes("multipart/form-data")) {
        return res.status(400).json({
          error: "Expected multipart/form-data containing a FASTA file"
        });
      }

      headers["Content-Type"] = contentType;
    }

    const options = {
      method: req.method,
      headers
    };

    if (req.method === "POST") {
      options.body = req;
      options.duplex = "half";
    }

    const upstream = await fetch(upstreamUrl, options);

    // Forward PLSDB response without changing its JSON structure.
    const responseBody = await upstream.arrayBuffer();

    res.status(upstream.status);

    res.setHeader(
      "Content-Type",
      upstream.headers.get("content-type") || "application/json"
    );

    return res.send(Buffer.from(responseBody));

  } catch (error) {
    console.error("PLSDB proxy error:", error);

    return res.status(502).json({
      error: "Could not communicate with the PLSDB API",
      details: error.message
    });
  }
}
