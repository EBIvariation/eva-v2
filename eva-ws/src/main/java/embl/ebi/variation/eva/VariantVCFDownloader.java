/*
 * Copyright 2016 EMBL - European Bioinformatics Institute
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package embl.ebi.variation.eva;


import embl.ebi.variation.eva.vcfDump.VariantExporterController;
import org.opencb.datastore.core.ObjectMap;
import org.opencb.datastore.core.QueryOptions;
import org.opencb.datastore.core.QueryResponse;
import org.opencb.datastore.core.QueryResult;
import org.opencb.opencga.lib.common.Config;
import org.opencb.opencga.storage.core.StorageManagerException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.servlet.http.HttpServletRequest;
import javax.ws.rs.*;
import javax.ws.rs.core.*;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Properties;

/**
 * Created by jmmut on 2016-01-22.
 *
 * @author Jose Miguel Mut Lopez &lt;jmmut@ebi.ac.uk&gt;
 */
@Path("/")
public class VariantVCFDownloader {

    private static final Logger logger = LoggerFactory.getLogger(VariantVCFDownloader.class);
    private final UriInfo uriInfo;
    private final MultivaluedMap<String, String> params;
    private final QueryOptions queryOptions;
    private final String sessionIp;

    protected String version;

    // Common output members
    protected long startTime;
    protected long endTime;
    private final String PROPERTIES = "/ws.properties";

//    protected static ObjectWriter jsonObjectWriter;
//    protected static ObjectMapper jsonObjectMapper;


//    {
//        jsonObjectMapper = new ObjectMapper();
//        jsonObjectWriter = jsonObjectMapper.writer();
//    }

    public VariantVCFDownloader(@PathParam("version") String version, @Context UriInfo uriInfo, @Context HttpServletRequest httpServletRequest) throws IOException {
        this.startTime = System.currentTimeMillis();
        this.version = version;
        this.uriInfo = uriInfo;
        this.params = uriInfo.getQueryParameters();
        logger.debug("requested: " + uriInfo.getRequestUri().toString());
        this.queryOptions = null;
        this.sessionIp = httpServletRequest.getRemoteAddr();
        Config.setOpenCGAHome(System.getenv("OPENCGA_HOME") != null ? System.getenv("OPENCGA_HOME") : "/opt/opencga");

    }

    @GET
    @Path("/echo/{message}")
    @Produces("text/plain")
    public Response echoGet(@PathParam("message") String message) {
        return buildResponse(Response.ok(message));
    }

    /**
     * example petition for downloading a vcf with some filters:
     * http://localhost:8080/eva-ws/rest/attach?species=hsapiens&assembly=grch37&studyIds=7&fileIds=5&maf=%3C0.3&region=20:60000:65000
     */
    @GET
    @Path("/attach")
    @Produces(MediaType.APPLICATION_OCTET_STREAM)
    public Response attach(@QueryParam("species") String species,
                           @QueryParam("assembly") String assembly,
                           @QueryParam("studyIds") String studyIds,
                           @QueryParam("fileIds") String fileIds)
            throws ClassNotFoundException, StorageManagerException, URISyntaxException, InstantiationException, IllegalAccessException, IOException {

        // TODO: parameter validation and show usage if error

//        startTime = System.currentTimeMillis();

        Properties properties = new Properties();
        InputStream propertiesStream = getClass().getResourceAsStream(PROPERTIES);
        if (propertiesStream == null) {
            throw new IOException("properties file " + PROPERTIES + " not found. VariantExporter can not continue.");
        }
        properties.load(propertiesStream);


        MultivaluedMap<String, String> queryParameters = uriInfo.getQueryParameters();
        logger.debug("queryParameters: " + queryParameters.toString());
        VariantExporterController vec = new VariantExporterController(
                species,
                getDbName(species, assembly),
                Arrays.asList(studyIds.split(",")),
                Arrays.asList(fileIds.split(",")),
                properties.getProperty("VCFDownloaderOutputFolder"),
                queryParameters);

        final List<String> files = vec.run();
        logger.info("returning dumped vcf as \"" + files.get(0) + "\".");

        return createOkResponse(new File(files.get(0)), MediaType.APPLICATION_OCTET_STREAM_TYPE, files.get(0));
    }

    /**
     * TODO: where should this go? This doesn't seem a good place
     */
    private String getDbName(String species, String assembly) {
        return "eva_" + species + "_" + assembly;
    }

    protected Response createJsonResponse(Object object) {
//        try {
        logger.info("on createJsonResponse, object=" + object.toString());
//            String entity = jsonObjectWriter.writeValueAsString(object);  // doesn't work, glassfish cannot instanciate VariantVCFDownloader
        Response.ResponseBuilder ok = Response.ok("{\"msg\":\"" + object.toString() + "\"}", MediaType.APPLICATION_JSON_TYPE);
        Response response = buildResponse(ok);
        return response;
//        } catch (JsonProcessingException e) {
//            System.out.println("object = " + object);
//            System.out.println("((QueryResponse)object).getResponse() = " + ((QueryResponse) object).getResponse());
//
//            System.out.println("e = " + e);
//            System.out.println("e.getMessage() = " + e.getMessage());
//            return createErrorResponse("Error parsing QueryResponse object:\n" + Arrays.toString(e.getStackTrace()));
//        }
    }

    protected Response createOkResponse(Object obj) {
        QueryResponse queryResponse = new QueryResponse();
        endTime = System.currentTimeMillis() - startTime;
        queryResponse.setTime(new Long(endTime - startTime).intValue());
        queryResponse.setApiVersion(version);

        // Guarantee that the QueryResponse object contains a list of results
        List list;
        if (obj instanceof List) {
            list = (List) obj;
        } else {
            list = new ArrayList();
            list.add(obj);
        }
        queryResponse.setResponse(list);

        return buildResponse(Response.ok());


    }

    protected Response createErrorResponse(Object o) {
        QueryResult<ObjectMap> result = new QueryResult();
        result.setErrorMsg(o.toString());
        return  createOkResponse(result);
    }

    protected Response createOkResponse(Object o1, MediaType o2) {
        return buildResponse(Response.ok(o1, o2));
    }

    protected Response createOkResponse(Object o1, MediaType o2, String fileName) {
        return buildResponse(Response.ok(o1, o2).header("content-disposition", "attachment; filename =" + fileName));
    }

    protected Response buildResponse(Response.ResponseBuilder responseBuilder) {
        endTime = System.currentTimeMillis();
        logger.debug("response completed in " + Long.toString(endTime - startTime) + " ms");
        return responseBuilder.header("Access-Control-Allow-Origin", "*").header("Access-Control-Allow-Headers", "x-requested-with, content-type").build();
    }
}
